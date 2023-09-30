use glam::{Vec3, DVec3};
use shalrath::parser::repr::parse_map;
use std::fs::File;
use unreal_asset::{
    base::types::PackageIndex,
    engine_version::EngineVersion,
    exports::{Export, ExportBaseTrait, ExportNormalTrait},
    properties::Property,
    Asset,
};

const GLOBAL_SCALE: f64 = 4.0;
const PI: f32 = std::f32::consts::PI;
const UNIT_X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const UNIT_Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const UNIT_Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

#[derive(Debug)]
struct UBox {
    position: DVec3,
    scale: DVec3,
    euler: DVec3,
}

impl UBox {
    fn new(p: DVec3, s: DVec3, r: DVec3) -> UBox {
        UBox {
            position: p,
            scale: s,
            euler: r,
        }
    }
}

fn main() {
    let mut args = std::env::args();
    _ = args.next();
    let input_path = args
        .next()
        .expect("Please specify map file as first argument.");
    let output_path = args
        .next()
        .expect("Please specify output directory as second argument.");
    let map_string = std::fs::read_to_string(input_path).expect("Failed to read map file.");
    let (_, map_ast) = parse_map(&map_string).expect("Failed to parse map.");
    let mut uboxes = vec![];
    for entity in map_ast.0 {
        for brush in entity.brushes.0 {
            let mut centroid = DVec3::new(0.0, 0.0, 0.0);
            let mut min = DVec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut max = DVec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for brush_plane in &brush.0 {
                let x = brush_plane.plane.v0.x as f64;
                let y = -brush_plane.plane.v0.y as f64;
                let z = brush_plane.plane.v0.z as f64;
                centroid.x += x;
                centroid.y += y;
                centroid.z += z;
                min.x = f64::min(x, min.x);
                min.y = f64::min(y, min.y);
                min.z = f64::min(z, min.z);
                max.x = f64::max(x, max.x);
                max.y = f64::max(y, max.y);
                max.z = f64::max(z, max.z);
            }
            let euler = derive_rotation(&brush) / PI * 180.0;
            uboxes.push(UBox::new(centroid / 6.0, max - min, DVec3::new(0.0, 0.0, 0.0)));
        }
    }

    main2(uboxes, &output_path);
}

fn main2(uboxes: Vec<UBox>, output_path: &str) {
    let data_file = std::io::BufReader::new(File::open("empty_level.umap").unwrap());
    let bulk_file = std::io::BufReader::new(File::open("empty_level.uexp").unwrap());

    let mut asset = Asset::new(data_file, Some(bulk_file), EngineVersion::VER_UE5_1, None).unwrap();

    for ubox in uboxes {
        add_box(&mut asset, ubox.position, ubox.scale, ubox.euler);
    }

    {
        let new_main_export_name = asset.add_fname("Slot1");
        let main_export = asset.get_export_mut(PackageIndex::new(6)).unwrap();
        assert_eq!(
            main_export
                .get_base_export_mut()
                .object_name
                .get_owned_content(),
            "empty_level"
        );
        main_export.get_base_export_mut().object_name = new_main_export_name;
    }

    if let Export::LevelExport(level_export) = asset.get_export_mut(PackageIndex::new(1)).unwrap() {
        for actor in &mut level_export.actors {
            if actor.index == 4 {
                actor.index = 0;
                println!("Removed original StaticMeshActor from PersistentLevel");
            }
        }
    } else {
        panic!("PersistentLevel not found");
    }

    let mut out_uasset = File::create(format!("{}/Slot1.umap", output_path)).unwrap();
    let mut out_uexp = File::create(format!("{}/Slot1.uexp", output_path)).unwrap();
    let _ = asset.write_data(&mut out_uasset, Some(&mut out_uexp));
}

fn add_box(asset: &mut Asset<std::io::BufReader<File>>, pos: DVec3, scale: DVec3, euler: DVec3) {
    let mut exp1 = asset.get_export(PackageIndex::new(4)).unwrap().clone();
    let exp2 = asset.get_export(PackageIndex::new(5)).unwrap().clone();
    asset.asset_data.exports.push(exp2);
    let idx2 = PackageIndex::new(asset.asset_data.exports.len() as i32);
    exp1.get_base_export_mut().object_name =
        asset.add_fname(&format!("StaticMeshActor_{}", idx2.index + 1));
    exp1.get_base_export_mut()
        .create_before_serialization_dependencies[0] = idx2;
    for prop in exp1.get_normal_export_mut().unwrap().properties.iter_mut() {
        if let Property::ObjectProperty(prop) = prop {
            prop.value = idx2;
        }
    }
    asset.asset_data.exports.push(exp1);
    let idx1 = PackageIndex::new(asset.asset_data.exports.len() as i32);
    let exp2 = asset.get_export_mut(idx2).unwrap();
    exp2.get_base_export_mut().outer_index = idx1;
    exp2.get_base_export_mut().create_before_create_dependencies[0] = idx1;
    let mut did_modify_scale = false;
    let mut did_modify_location = false;
    let mut did_modify_rotation = false;
    for prop in exp2.get_normal_export_mut().unwrap().properties.iter_mut() {
        if let Property::StructProperty(prop) = prop {
            if prop.name.get_owned_content() == "RelativeLocation" {
                for prop in &mut prop.value {
                    if let Property::VectorProperty(prop) = prop {
                        assert!(!did_modify_location);
                        prop.value.x.0 = GLOBAL_SCALE * pos.x;
                        prop.value.y.0 = GLOBAL_SCALE * pos.y;
                        prop.value.z.0 = GLOBAL_SCALE * pos.z;
                        did_modify_location = true;
                    }
                }
            }
            if prop.name.get_owned_content() == "RelativeScale3D" {
                for prop in &mut prop.value {
                    if let Property::VectorProperty(prop) = prop {
                        assert!(!did_modify_scale);
                        prop.value.x.0 = GLOBAL_SCALE * scale.x;
                        prop.value.y.0 = GLOBAL_SCALE * scale.y;
                        prop.value.z.0 = GLOBAL_SCALE * scale.z;
                        did_modify_scale = true;
                    }
                }
            }
            if prop.name.get_owned_content() == "RelativeRotation" {
                for prop in &mut prop.value {
                    if let Property::RotatorProperty(prop) = prop {
                        assert!(!did_modify_rotation);
                        prop.value.x.0 = euler.y;
                        prop.value.y.0 = euler.x;
                        prop.value.z.0 = euler.z;
                        did_modify_rotation = true;
                    }
                }
            }
        }
    }
    assert!(did_modify_scale);
    assert!(did_modify_location);
    assert!(did_modify_rotation);

    let level_export = asset.get_export_mut(PackageIndex::new(1)).unwrap();
    level_export
        .get_base_export_mut()
        .create_before_serialization_dependencies
        .push(idx1);
    if let Export::LevelExport(level_export) = level_export {
        level_export.actors.push(idx1);
    }
}

fn vec3(p: shalrath::repr::Point) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}
fn dvec3(v: Vec3) -> DVec3 {
    DVec3::new(v.x as f64, v.y as f64, v.z as f64)
}
fn calculate_normal(plane: shalrath::repr::TrianglePlane) -> Vec3 {
    let p0 = vec3(plane.v0);
    let p1 = vec3(plane.v1);
    let p2 = vec3(plane.v2);
    let a = p2 - p0;
    let b = p1 - p0;
    a.cross(b)
}
fn angle_between(a: &Vec3, b: &Vec3) -> f32 {
    trunc_with_precision(a.dot(*b) / (a.length() * b.length()), 1e-5).acos()
}
fn trunc_with_precision(f: f32, precision: f32) -> f32 {
    (f / precision).trunc() * precision
}
fn approx_eq(a: f32, b: f32, epsilon: f32) -> bool {
    (a - b).abs() <= epsilon
}
fn is_box(brush: &shalrath::repr::Brush) -> bool {
    if brush.0.len() != 6 {
        return false;
    }
    let normals: Vec<Vec3> = brush.0.iter().map(|plane| calculate_normal(plane.plane)).collect();
    // TODO this could be more efficient
    for a in &normals {
        let mut angle0 = 0;
        let mut angle90 = 0;
        let mut angle180 = 0;
        for b in &normals {
            let angle = angle_between(a, b).abs();
            if approx_eq(angle, 0.0, 1e-2) {
                angle0 += 1;
            } else if approx_eq(angle, 0.5*PI, 1e-2) {
                angle90 += 1;
            } else if approx_eq(angle, PI, 1e-2) {
                angle180 += 1;
            }
        }
        if angle0 != 1 || angle90 != 4 || angle180 != 1 {
            return false;
        }
    }
    true
}
fn derive_rotation(brush: &shalrath::repr::Brush) -> Vec3 {
    assert!(is_box(brush));
    let normals: Vec<Vec3> = brush.0.iter().map(|plane| calculate_normal(plane.plane)).collect();
    let mut is_aligned_x = false;
    let mut is_aligned_y = false;
    let mut is_aligned_z = false;
    let mut is_aligned_normal = [false; 6];
    for (i, n) in normals.iter().enumerate() {
        if n.x != 0.0 && n.y == 0.0 && n.z == 0.0 {
            is_aligned_x = true;
            is_aligned_normal[i] = true;
        }
        if n.x == 0.0 && n.y != 0.0 && n.z == 0.0 {
            is_aligned_y = true;
            is_aligned_normal[i] = true;
        }
        if n.x == 0.0 && n.y == 0.0 && n.z != 0.0 {
            is_aligned_z = true;
            is_aligned_normal[i] = true;
        }
    }
    if is_aligned_x && is_aligned_y && is_aligned_z {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    if is_aligned_x || is_aligned_y || is_aligned_z {
        let mut result = Vec3::new(0.0, 0.0, 0.0);
        let (ref_axis, dst) = if is_aligned_x {
            println!("x aligned");
            (&UNIT_Z, &mut result.x)
        } else if is_aligned_y {
            println!("y aligned");
            (&UNIT_X, &mut result.y)
        } else {
            println!("z aligned");
            (&UNIT_Y, &mut result.z)
        };
        let mut pos_quadrant_normal_found = false;
        for (i, n) in normals.iter().enumerate() {
            if is_aligned_normal[i] {
                continue;
            }
            if is_in_positive_quadrant(n) {
                pos_quadrant_normal_found = true;
                
                *dst = angle_between(n, ref_axis);
            }
        }
        assert!(pos_quadrant_normal_found);
        return result;
    }
    // need to rotate on all 3 axes
    panic!("can't handle 3-axis euler rotations right now :(. please stick to single axis rotated boxes.");
    //Vec3::new(0.0, 0.0, 0.0)
}
fn is_in_positive_quadrant(v: &Vec3) -> bool {
    v.x >= 0.0 && v.y >= 0.0 && v.z >= 0.0
}
