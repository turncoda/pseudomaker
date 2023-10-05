use clap::Parser;
use glam::{DQuat, DVec3, Mat3, Vec3};
use itertools::Itertools;
use shalrath::parser::repr::parse_map;
use std::fs::File;
use std::io::Cursor;
use unreal_asset::{
    base::types::PackageIndex,
    engine_version::EngineVersion,
    exports::{Export, ExportBaseTrait, ExportNormalTrait},
    properties::Property,
    Asset,
};

/// Convert simple Quake maps into Pseudoregalia levels.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to Quake map file.
    #[arg(short, long)]
    map: String,

    /// Asset path in pak. (Example: "/Game/Maps/MyMap")
    #[arg(short, long)]
    path: String,

    /// Path of output pak file.
    #[arg(short, long)]
    out: String,

    /// Scale of level.
    #[arg(short, long, default_value_t = 1.0)]
    scale: f64,
}

// make this a command line arg or world settings param
const GLOBAL_SCALE: f64 = 4.0;
const PI: f32 = std::f32::consts::PI;
const PI64: f64 = std::f64::consts::PI;
const UNIT_X: DVec3 = DVec3::new(1.0, 0.0, 0.0);
const UNIT_Y: DVec3 = DVec3::new(0.0, 1.0, 0.0);
const UNIT_Z: DVec3 = DVec3::new(0.0, 0.0, 1.0);
const U1_CUBE_UASSET: &[u8] = include_bytes!("u1_cube.uasset");
const U1_CUBE_UBULK: &[u8] = include_bytes!("u1_cube.ubulk");
const U1_CUBE_UEXP: &[u8] = include_bytes!("u1_cube.uexp");
const EMPTY_LEVEL_UMAP: &[u8] = include_bytes!("empty_level.umap");
const EMPTY_LEVEL_UEXP: &[u8] = include_bytes!("empty_level.uexp");

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
    let args = Args::parse();
    let relative_path = args
        .path
        .strip_prefix("/Game/")
        .expect("Please pass in a path starting with \"/Game/\".");
    let export_name = args.path.split("/").last().unwrap();
    let map_string = std::fs::read_to_string(args.map).expect("Failed to read map file.");
    let (_, map_ast) = parse_map(&map_string).expect("Failed to parse map.");
    let mut uboxes = vec![];
    for entity in map_ast.0 {
        for brush in entity.brushes.0 {
            // TrenchBroom uses a right-handed coordinate system with +z facing up.
            // Unreal uses a left-handed coordinate system with +z facing up.
            // So, we mirror across the xz-plane when converting.
            let brush = mirror_xz(&brush);
            let vertices = calculate_vertices(&brush);
            let centroid = dvec3(avg(&vertices));
            let rot = derive_rotation(&brush);
            let mut min = DVec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut max = DVec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for p in vertices {
                let p = dvec3(p);
                let rp = rot.inverse() * (p - centroid);
                let x = rp.x as f64;
                let y = rp.y as f64;
                let z = rp.z as f64;
                min.x = f64::min(x, min.x);
                min.y = f64::min(y, min.y);
                min.z = f64::min(z, min.z);
                max.x = f64::max(x, max.x);
                max.y = f64::max(y, max.y);
                max.z = f64::max(z, max.z);
            }
            let (z, y, x) = rot.to_euler(glam::EulerRot::ZYX);
            let mut euler = DVec3::new(-y, z, -x);
            euler = euler / PI64 * 180.0;
            uboxes.push(UBox::new(centroid, max - min, euler));
        }
    }

    let mut asset = Asset::new(
        Cursor::new(EMPTY_LEVEL_UMAP),
        Some(Cursor::new(EMPTY_LEVEL_UEXP)),
        EngineVersion::VER_UE5_1,
        None,
    )
    .unwrap();

    for ubox in uboxes {
        add_box(
            &mut asset,
            ubox.position,
            ubox.scale,
            ubox.euler,
            args.scale,
        );
    }

    {
        let new_main_export_name = asset.add_fname(&export_name);
        // TODO don't hardcode these export indices
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
            }
        }
    } else {
        panic!("PersistentLevel not found");
    }

    let mut umap_data = vec![];
    let mut uexp_data = vec![];
    let _ = asset.write_data(
        &mut Cursor::new(&mut umap_data),
        Some(&mut Cursor::new(&mut uexp_data)),
    );

    let mut pak = repak::PakWriter::new(
        File::create(args.out).unwrap(),
        repak::Version::V11,
        "../../../pseudoregalia/Content/".to_string(),
        None,
    );
    pak.write_file("Mods/PseudoMaker/u1_cube.uasset", U1_CUBE_UASSET)
        .unwrap();
    pak.write_file("Mods/PseudoMaker/u1_cube.ubulk", U1_CUBE_UBULK)
        .unwrap();
    pak.write_file("Mods/PseudoMaker/u1_cube.uexp", U1_CUBE_UEXP)
        .unwrap();
    pak.write_file(&format!("{}.umap", relative_path), &umap_data)
        .unwrap();
    pak.write_file(&format!("{}.uexp", relative_path), &uexp_data)
        .unwrap();
    pak.write_index().unwrap();
}

fn add_box(
    asset: &mut Asset<std::io::Cursor<&[u8]>>,
    pos: DVec3,
    scale: DVec3,
    euler: DVec3,
    scale_arg: f64,
) {
    // TODO don't hardcode these export indices
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
                        prop.value.x.0 = scale_arg * GLOBAL_SCALE * pos.x;
                        prop.value.y.0 = scale_arg * GLOBAL_SCALE * pos.y;
                        prop.value.z.0 = scale_arg * GLOBAL_SCALE * pos.z;
                        did_modify_location = true;
                    }
                }
            }
            if prop.name.get_owned_content() == "RelativeScale3D" {
                for prop in &mut prop.value {
                    if let Property::VectorProperty(prop) = prop {
                        assert!(!did_modify_scale);
                        prop.value.x.0 = scale_arg * GLOBAL_SCALE * scale.x;
                        prop.value.y.0 = scale_arg * GLOBAL_SCALE * scale.y;
                        prop.value.z.0 = scale_arg * GLOBAL_SCALE * scale.z;
                        did_modify_scale = true;
                    }
                }
            }
            if prop.name.get_owned_content() == "RelativeRotation" {
                for prop in &mut prop.value {
                    if let Property::RotatorProperty(prop) = prop {
                        assert!(!did_modify_rotation);
                        prop.value.x.0 = euler.x;
                        prop.value.y.0 = euler.y;
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
    let normals: Vec<Vec3> = brush
        .0
        .iter()
        .map(|plane| calculate_normal(plane.plane))
        .collect();
    // TODO this could be more efficient
    for a in &normals {
        let mut angle0 = 0;
        let mut angle90 = 0;
        let mut angle180 = 0;
        for b in &normals {
            let angle = angle_between(a, b).abs();
            if approx_eq(angle, 0.0, 1e-2) {
                angle0 += 1;
            } else if approx_eq(angle, 0.5 * PI, 1e-2) {
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
fn avg(pts: &[Vec3]) -> Vec3 {
    let sum = pts.iter().fold(Vec3::new(0.0, 0.0, 0.0), |a, b| a + *b);
    sum / pts.len() as f32
}
fn calculate_vertices(brush: &shalrath::repr::Brush) -> Vec<Vec3> {
    let mut plane_intersections = vec![];
    let planes_and_normals: Vec<_> = brush
        .0
        .iter()
        .map(|plane| (plane.plane, calculate_normal(plane.plane)))
        .collect();
    for c in planes_and_normals.iter().combinations(3) {
        let (p0, n0) = c[0];
        let (p1, n1) = c[1];
        let (p2, n2) = c[2];
        let a0 = angle_between(n0, n1).abs();
        let a1 = angle_between(n1, n2).abs();
        let a2 = angle_between(n2, n0).abs();
        let b0 = approx_eq(a0, 0.5 * PI, 1e-2);
        let b1 = approx_eq(a1, 0.5 * PI, 1e-2);
        let b2 = approx_eq(a2, 0.5 * PI, 1e-2);
        if b0 && b1 && b2 {
            let x0 = vec3(p0.v0);
            let x1 = vec3(p1.v0);
            let x2 = vec3(p2.v0);
            let u0 = n0.normalize();
            let u1 = n1.normalize();
            let u2 = n2.normalize();
            let det = Mat3::from_cols(u0, u1, u2).determinant();
            let c0 = x0.dot(u0) * u1.cross(u2);
            let c1 = x1.dot(u1) * u2.cross(u0);
            let c2 = x2.dot(u2) * u0.cross(u1);
            let intersection = det.powi(-1) * (c0 + c1 + c2);
            plane_intersections.push(intersection);
        }
    }
    assert_eq!(8, plane_intersections.len());
    plane_intersections
}
fn point(v: DVec3) -> shalrath::repr::Point {
    shalrath::repr::Point {
        x: v.x as f32,
        y: v.y as f32,
        z: v.z as f32,
    }
}
fn rotate_brush(brush: &shalrath::repr::Brush, rot: DQuat) -> shalrath::repr::Brush {
    let centroid = avg(&calculate_vertices(&brush));
    let mut planes = vec![];
    for plane in &brush.0 {
        let mut plane = plane.clone();
        plane.plane.v0 = point(rot * dvec3(vec3(plane.plane.v0) - centroid) + dvec3(centroid));
        plane.plane.v1 = point(rot * dvec3(vec3(plane.plane.v1) - centroid) + dvec3(centroid));
        plane.plane.v2 = point(rot * dvec3(vec3(plane.plane.v2) - centroid) + dvec3(centroid));
        planes.push(plane);
    }
    shalrath::repr::Brush(planes)
}
fn derive_rotation(brush: &shalrath::repr::Brush) -> DQuat {
    assert!(is_box(brush));
    let normals: Vec<Vec3> = brush
        .0
        .iter()
        .map(|plane| calculate_normal(plane.plane))
        .collect();
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
        return DQuat::IDENTITY;
    }
    let unaligned_normals: Vec<_> = normals
        .iter()
        .enumerate()
        .filter_map(|(i, n)| {
            if is_aligned_normal[i] {
                None
            } else {
                Some(n.clone())
            }
        })
        .collect();
    if unaligned_normals.len() == 4 {
        let pqn = find_pos_quadrant_normal(&unaligned_normals)
            .unwrap()
            .normalize();
        // using right-handed coordinate system
        if is_aligned_x {
            return DQuat::from_rotation_arc(UNIT_Y, dvec3(pqn));
        }
        if is_aligned_y {
            return DQuat::from_rotation_arc(UNIT_Z, dvec3(pqn));
        }
        if is_aligned_z {
            return DQuat::from_rotation_arc(UNIT_X, dvec3(pqn));
        }
        panic!("unreachable");
    }
    panic!("can't handle multi-axis rotations right now. sorry :/");
    let phn = find_pos_hemisphere_normal(&unaligned_normals)
        .unwrap()
        .normalize();
    let z_aligning_rotation = DQuat::from_rotation_arc(UNIT_Z, dvec3(phn));
    z_aligning_rotation + derive_rotation(&rotate_brush(&brush, z_aligning_rotation.inverse()))
}
fn find_pos_hemisphere_normal(normals: &[Vec3]) -> Option<Vec3> {
    for n in normals {
        if n.z >= 0.0 {
            return Some(n.clone());
        }
    }
    None
}
fn find_pos_quadrant_normal(normals: &[Vec3]) -> Option<Vec3> {
    for n in normals {
        if is_in_positive_quadrant(n) {
            return Some(n.clone());
        }
    }
    None
}
fn is_in_positive_quadrant(v: &Vec3) -> bool {
    let epsilon = 1e-3;
    v.x >= -epsilon && v.y >= -epsilon && v.z >= -epsilon
}
fn flip_y(point: &mut shalrath::repr::Point) {
    point.y = -point.y;
}
fn mirror_xz(brush: &shalrath::repr::Brush) -> shalrath::repr::Brush {
    let mut planes = vec![];
    for plane in &brush.0 {
        let mut mirrored = plane.clone();
        flip_y(&mut mirrored.plane.v0);
        flip_y(&mut mirrored.plane.v1);
        flip_y(&mut mirrored.plane.v2);
        planes.push(mirrored);
    }
    shalrath::repr::Brush::new(planes)
}
