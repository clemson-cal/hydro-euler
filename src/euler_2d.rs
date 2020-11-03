use crate::geometry::{Direction, Vector3d};


// ============================================================================
#[derive(Copy, Clone, Debug)]
pub struct Conserved(pub f64, pub f64, pub f64, pub f64);

#[derive(Copy, Clone, Debug)]
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
impl Into<[f64; 4]> for Primitive { fn into(self) -> [f64; 4] { [self.0, self.1, self.2, self.3] } }
impl From<[f64; 4]> for Primitive { fn from(a:  [f64; 4]) -> Primitive { Primitive(a[0], a[1], a[2], a[3]) } }

impl Into<[f64; 4]> for Conserved { fn into(self) -> [f64; 4] { [self.0, self.1, self.2, self.3] } }
impl From<[f64; 4]> for Conserved { fn from(a:  [f64; 4]) -> Conserved { Conserved(a[0], a[1], a[2], a[3]) } }

impl Default for Conserved { fn default() -> Self { Conserved(0.0, 0.0, 0.0, 0.0) } }
impl Default for Primitive { fn default() -> Self { Primitive(0.0, 0.0, 0.0, 0.0) } }

impl Conserved { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e} }
impl Primitive { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e && self.3.abs() < e} }




// ============================================================================
impl Conserved {
    pub fn mass_density     (self)  -> f64 { self.0 }
    pub fn momentum_1       (self)  -> f64 { self.1 }
    pub fn momentum_2       (self)  -> f64 { self.2 }
    pub fn energy_density   (self)  -> f64 { self.3 }

    pub fn momentum_vector  (self)  -> Vector3d {
        Vector3d(self.momentum_1(), self.momentum_2(), 0.0)
    }
    pub fn momentum(self, direction: Direction) -> f64 {
        match direction{
            Direction::X => self.momentum_1(),
            Direction::Y => self.momentum_2(),
            Direction::Z => 0.0,
        }
    }
    pub fn momentum_squared(self) -> f64 {
        let s1 = self.momentum_1();
        let s2 = self.momentum_2();
        s1*s1 + s2*s2 //Note: We use the normalized coordinate basis, so this is the correct way to square a vector in any coordinates.
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
    pub fn mass_density(self) -> f64 { self.0 }
    pub fn velocity_1  (self) -> f64 { self.1 }
    pub fn velocity_2  (self) -> f64 { self.2 }
    pub fn gas_pressure(self) -> f64 { self.3 }

    pub fn velocity    (self, direction: Direction) -> f64 {
        match direction {
            Direction::X => self.velocity_1(),
            Direction::Y => self.velocity_2(),
            Direction::Z => 0.0,
        }
    }

    pub fn sound_speed_squared(self, gamma_law_index: f64) -> f64 {
        gamma_law_index * self.gas_pressure() / self.mass_density()
    }

    pub fn outer_wavespeeds(self, direction: Direction, gamma_law_index: f64) -> (f64, f64) {
        let cs = self.sound_speed_squared(gamma_law_index).sqrt();
        let vn = self.velocity(direction);
        (vn - cs, vn + cs)
    }

//    pub fn max_wavespeed(self, gamma_law_index: f64) -> f64 {
//        //let (am, ap) = self.outer_wavespeeds(gamma_law_index);
//        let (alm, alp) = pl.outer_wavespeeds(direction, gamma_law_index);
//        let (arm, arp) = pr.outer_wavespeeds(direction, gamma_law_index);
//        f64::max(am.abs(), ap.abs())
//    }

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

    pub fn reflect(self) -> Primitive {
        Primitive(self.0, -self.1, -self.2, self.3)
    }
}




// ============================================================================
impl Primitive
{
    pub fn spherical_geometry_source_terms(self, spherical_radius: f64, polar_angle_theta: f64) -> Conserved
    {
        let v_theta       = &self.velocity_2;
        let v_phi         = 0.0;
        let pressure_term = 2.0 * self.gas_pressure() / spherical_radius;
        let velocity_term = self.density() * (v_theta * v_theta + v_phi * v_phi) / spherical_radius;
        let theta_term    = pressure_term + velocity_term
        Conserved(0.0, theta_term, 0.0, 0.0)
    }
}




//=============================================================================
pub enum RiemannSolverMode
{
    HlleFlux,
    HlleFluxAcrossMovingFace(f64)
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




// ============================================================================
#[cfg(test)]
mod tests
{
	use crate::euler_1d::Primitive;

    fn panic_unless_recovery_is_accurate(primitive: Primitive)
    {
        let gamma_law_index = 5.0 / 3.0;
        let u = primitive.to_conserved(gamma_law_index);
        let p = u.to_primitive(gamma_law_index).unwrap();
        println!("{:?} {:?}", primitive, p);
        assert!(f64::abs(primitive.velocity_1() - p.velocity_1()) < 1e-10);
        assert!(f64::abs(primitive.velocity_2() - p.velocity_2()) < 1e-10);
        assert!(f64::abs(primitive.gas_pressure() - p.gas_pressure()) / primitive.gas_pressure() < 1e-10);
    }

    #[test]
    fn can_recover_primitive()
    {
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.2, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 0.2, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.5, 0.5, 1e-3));
        panic_unless_recovery_is_accurate(Primitive(1.0, 5.0, 5.0, 1e+3));
    }
}
