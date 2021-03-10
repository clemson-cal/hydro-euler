use crate::geometry::{Direction, Vector3d};




// ============================================================================
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
#[cfg_attr(feature="serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Conserved(pub f64, pub f64, pub f64, pub f64);

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
#[cfg_attr(feature="serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Primitive(pub f64, pub f64, pub f64, pub f64);




// ============================================================================
impl std::ops::Add<Primitive> for Primitive { type Output = Self; fn add(self, u: Primitive) -> Primitive { Primitive(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl std::ops::Sub<Primitive> for Primitive { type Output = Self; fn sub(self, u: Primitive) -> Primitive { Primitive(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl std::ops::Mul<f64> for Primitive { type Output = Primitive; fn mul(self, a: f64) -> Primitive { Primitive(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl std::ops::Div<f64> for Primitive { type Output = Primitive; fn div(self, a: f64) -> Primitive { Primitive(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a, self.3 / a) } }




// ============================================================================
impl From<Primitive> for [f64; 4] { fn from(a: Primitive) -> [f64; 4] { [a.0, a.1, a.2, a.3] } }
impl From<Conserved> for [f64; 4] { fn from(a: Conserved) -> [f64; 4] { [a.0, a.1, a.2, a.3] } }
impl From<[f64; 4]> for Primitive { fn from(a: [f64; 4]) -> Primitive { Primitive(a[0], a[1], a[2], a[3]) } }
impl From<[f64; 4]> for Conserved { fn from(a: [f64; 4]) -> Conserved { Conserved(a[0], a[1], a[2], a[3]) } }

impl Default for Conserved { fn default() -> Self { Conserved(0.0, 0.0, 0.0, 0.0) } }
impl Default for Primitive { fn default() -> Self { Primitive(0.0, 0.0, 0.0, 0.0) } }

impl Conserved { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e} }
impl Primitive { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e} }




// ============================================================================
impl Conserved {

    pub fn mass_density(self) -> f64 {
        self.0
    }

    pub fn momentum_1(self) -> f64 {
        self.1
    }

    pub fn momentum_2(self) -> f64 {
        self.2
    }

    pub fn energy_density(self) -> f64 {
        self.3
    }

    pub fn momentum_vector(self)  -> Vector3d {
        Vector3d(self.momentum_1(), self.momentum_2(), 0.0)
    }

    pub fn momentum(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.momentum_1(),
            Direction::Y => self.momentum_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn momentum_squared(self) -> f64 {
        let s1 = self.momentum_1();
        let s2 = self.momentum_2();
        s1 * s1 + s2 * s2
    }

    pub fn to_primitive(self, gamma_law_index: f64) -> Primitive {
    	let ek = 0.5 * self.momentum_squared() / self.mass_density();
    	let et = self.energy_density() - ek;
    	let pg = et * (gamma_law_index - 1.0);
    	let vx = self.momentum_1() / self.mass_density();
        let vy = self.momentum_2() / self.mass_density();
    	Primitive(self.mass_density(), vx, vy, pg)
    }
}




// ============================================================================
impl Primitive {

    pub fn mass_density(self) -> f64 {
        self.0
    }

    pub fn velocity_1(self) -> f64 {
        self.1
    }

    pub fn velocity_2(self) -> f64 {
        self.2
    }

    pub fn gas_pressure(self) -> f64 {
        self.3
    }

    pub fn velocity(self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.velocity_1(),
            Direction::Y => self.velocity_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn velocity_squared(self) -> f64 {
        self.1 * self.1 + self.2 * self.2
    }

    pub fn sound_speed_squared(self, gamma_law_index: f64) -> f64 {
        gamma_law_index * self.gas_pressure() / self.mass_density()
    }

    pub fn specific_kinetic_energy(self) -> f64 {
        0.5 * self.velocity_squared()
    }

    pub fn specific_internal_energy(self, gamma_law_index: f64) -> f64 {
        self.gas_pressure() / self.mass_density() / (gamma_law_index - 1.0)
    }

    pub fn mach_number(self, gamma_law_index: f64) -> f64 {
        (self.velocity_squared() / self.sound_speed_squared(gamma_law_index)).sqrt()
    }

    pub fn outer_wavespeeds(self, direction: Direction, gamma_law_index: f64) -> (f64, f64) {
        let cs = self.sound_speed_squared(gamma_law_index).sqrt();
        let vn = self.velocity(direction);
        (vn - cs, vn + cs)
    }

    pub fn max_signal_speed(self, gamma_law_index: f64) -> f64 {
        f64::sqrt(self.velocity_squared()) + f64::sqrt(self.sound_speed_squared(gamma_law_index))
    }

    pub fn to_conserved(self, gamma_law_index: f64) -> Conserved {
        let m   = self.mass_density();
        let p   = self.gas_pressure();
        let vsq = self.velocity_1().powi(2) + self.velocity_2().powi(2);

        Conserved(
            m,
            m * self.velocity_1(),
            m * self.velocity_2(),
            m * vsq * 0.5 + p / (gamma_law_index - 1.0)
        )
    }

    pub fn flux_vector(self, direction: Direction, gamma_law_index: f64) -> Conserved {
        let pg = self.gas_pressure();
        let vn = self.velocity(direction);

        let pressure_term = Conserved(
            0.0,
            pg * direction.along(Direction::X),
            pg * direction.along(Direction::Y),
            pg * vn);

        let advective_term = self.to_conserved(gamma_law_index) * vn;

        return advective_term + pressure_term;
    }

    pub fn reflect(self, direction: Direction) -> Primitive {
        match direction {
            Direction::X => Primitive(self.0, -self.1, self.2, self.3),
            Direction::Y => Primitive(self.0, self.1, -self.2, self.3),
            Direction::Z => panic!("no reflection for Direction::Z"),
        }
    }
}




// ============================================================================
impl Primitive {

    pub fn spherical_geometry_source_terms(self, spherical_radius: f64, polar_angle_theta: f64) -> Conserved {
        let cotq = f64::tan(std::f64::consts::FRAC_PI_2 - polar_angle_theta);
        let ur = self.velocity_1();
        let uq = self.velocity_2();
        let up = 0.0;
        let pg = self.gas_pressure();
        let d0 = self.mass_density();
        let sd = 0.0;
        let sr = (2.0  * pg + d0 * (uq * uq        + up * up)) / spherical_radius;
        let sq = (cotq * pg + d0 * (up * up * cotq - ur * uq)) / spherical_radius;
        // let sp =        -up * d0 * (ur + uq * cotq) / spherical_radius;
        let se = 0.0;
        Conserved(sd, sr, sq, se)
    }
}




// ============================================================================
pub fn riemann_hlle(pl: Primitive, pr: Primitive, direction: Direction, gamma_law_index: f64) -> Conserved {
    let ul = pl.to_conserved(gamma_law_index);
    let ur = pr.to_conserved(gamma_law_index);
    let fl = pl.flux_vector(direction, gamma_law_index);
    let fr = pr.flux_vector(direction, gamma_law_index);

    let (alm, alp) = pl.outer_wavespeeds(direction, gamma_law_index);
    let (arm, arp) = pr.outer_wavespeeds(direction, gamma_law_index);
    let ap = alp.max(arp).max(0.0);
    let am = alm.min(arm).min(0.0);

    (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am)
}





/**
 * Return Godunov fluxes of the conserved quantities and a passive scalar. The
 * inputs sl and sr are primitive-like: they are the volumetric scalar density,
 * just like the primitives contain the volumetric mass density.
 */
pub fn riemann_hlle_scalar(
    pl: Primitive,
    pr: Primitive,
    sl: f64,
    sr: f64,
    nhat: Direction,
    gamma_law_index: f64) -> (Conserved, f64)
{
    let cl = sl / pl.mass_density(); // scalar concentration
    let cr = sr / pr.mass_density();
    let f = riemann_hlle(pl, pr, nhat, gamma_law_index);

    let g = if f.mass_density() < 0.0 {
        cr * f.mass_density()
    } else {
        cl * f.mass_density()
    };
    (f, g)
}




// ============================================================================
#[cfg(test)]
mod tests {
	use crate::euler_2d::Primitive;

    fn panic_unless_recovery_is_accurate(primitive: Primitive) {
        let gamma_law_index = 5.0 / 3.0;
        let u = primitive.to_conserved(gamma_law_index);
        let p = u.to_primitive(gamma_law_index);
        println!("{:?} {:?}", primitive, p);
        assert!(f64::abs(primitive.velocity_1() - p.velocity_1()) < 1e-10);
        assert!(f64::abs(primitive.velocity_2() - p.velocity_2()) < 1e-10);
        assert!(f64::abs(primitive.gas_pressure() - p.gas_pressure()) / primitive.gas_pressure() < 1e-10);
    }

    #[test]
    fn can_recover_primitive() {
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.2, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.2, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.5, 0.5, 1e-3));
        panic_unless_recovery_is_accurate(Primitive(1.0, 5.0, 5.0, 1e+3));
    }
}
