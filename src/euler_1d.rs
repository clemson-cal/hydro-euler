// ============================================================================
#[derive(Clone, Copy, Debug)] #[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
pub struct Conserved(pub f64, pub f64, pub f64);

#[derive(Clone, Copy, Debug)] #[cfg_attr(feature="hdf5", repr(C), derive(hdf5::H5Type))]
pub struct Primitive(pub f64, pub f64, pub f64);




// ============================================================================
impl std::ops::Add<Primitive> for Primitive { type Output = Self; fn add(self, u: Primitive) -> Primitive { Primitive(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl std::ops::Sub<Primitive> for Primitive { type Output = Self; fn sub(self, u: Primitive) -> Primitive { Primitive(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl std::ops::Mul<f64> for Primitive { type Output = Primitive; fn mul(self, a: f64) -> Primitive { Primitive(self.0 * a, self.1 * a, self.2 * a) } }
impl std::ops::Div<f64> for Primitive { type Output = Primitive; fn div(self, a: f64) -> Primitive { Primitive(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
impl From<Primitive> for [f64; 3] { fn from(a: Primitive) -> [f64; 3] { [a.0, a.1, a.2] } }
impl From<Conserved> for [f64; 3] { fn from(a: Conserved) -> [f64; 3] { [a.0, a.1, a.2] } }
impl From<[f64; 3]> for Primitive { fn from(a: [f64; 3]) -> Primitive { Primitive(a[0], a[1], a[2]) } }
impl From<[f64; 3]> for Conserved { fn from(a: [f64; 3]) -> Conserved { Conserved(a[0], a[1], a[2]) } }

impl Default for Conserved { fn default() -> Self { Conserved(0.0, 0.0, 0.0) } }
impl Default for Primitive { fn default() -> Self { Primitive(0.0, 0.0, 0.0) } }

impl Conserved { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e } }
impl Primitive { pub fn small(self, e: f64) -> bool { self.0.abs() < e && self.1.abs() < e && self.2.abs() < e } }




// ============================================================================
impl Conserved {
    pub fn mass_density     (self)  -> f64 { self.0 }
    pub fn momentum_density (self)  -> f64 { self.1 }
    pub fn energy_density   (self)  -> f64 { self.2 }

    pub fn momentum_squared(self) -> f64 {
        self.momentum_density().powi(2)
    }

    pub fn to_primitive(self, gamma_law_index: f64) -> Primitive {
    	let ek = 0.5 * self.momentum_squared() / self.mass_density();
    	let et = self.energy_density() - ek;
    	let pg = et * (gamma_law_index - 1.0);
    	let vx = self.momentum_density() / self.mass_density();
    	Primitive(self.mass_density(), vx, pg)
    }
}




// ============================================================================
impl Primitive {
    pub fn mass_density(self) -> f64 { self.0 }
    pub fn velocity    (self) -> f64 { self.1 }
    pub fn gas_pressure(self) -> f64 { self.2 }

    pub fn sound_speed_squared(self, gamma_law_index: f64) -> f64 {
        gamma_law_index * self.gas_pressure() / self.mass_density()
    }

    pub fn outer_wavespeeds(self, gamma_law_index: f64) -> (f64, f64) {
        let cs = self.sound_speed_squared(gamma_law_index).sqrt();
        let vx = self.velocity();
        (vx - cs, vx + cs)
    }

    pub fn max_wavespeed(self, gamma_law_index: f64) -> f64 {
        let (am, ap) = self.outer_wavespeeds(gamma_law_index);
        f64::max(am.abs(), ap.abs())
    }

    pub fn to_conserved(self, gamma_law_index: f64) -> Conserved {
        let m = self.mass_density();
        let p = self.gas_pressure();

        Conserved(
            m,
            m * self.velocity(),
            m * self.velocity().powi(2) * 0.5 + p / (gamma_law_index - 1.0)
        )
    }

    pub fn flux_vector(self, gamma_law_index: f64) -> Conserved {
        let pg = self.gas_pressure();
        let vx = self.velocity();

        let pressure_term = Conserved(
            0.0,
            pg,
            pg * vx);

        let advective_term = self.to_conserved(gamma_law_index) * vx;

        return advective_term + pressure_term;
    }

    pub fn reflect(self) -> Primitive {
        Primitive(self.0, -self.1, self.2)
    }
}




// ============================================================================
impl Primitive
{
    pub fn spherical_geometry_source_terms(self, spherical_radius: f64) -> Conserved
    {
        Conserved(0.0, 2.0 * self.gas_pressure() / spherical_radius, 0.0)
    }
}




// ============================================================================
pub fn riemann_hlle(pl: Primitive, pr: Primitive, gamma_law_index: f64) -> Conserved {
    let ul = pl.to_conserved(gamma_law_index);
    let ur = pr.to_conserved(gamma_law_index);
    let fl = pl.flux_vector(gamma_law_index);
    let fr = pr.flux_vector(gamma_law_index);

    let (alm, alp) = pl.outer_wavespeeds(gamma_law_index);
    let (arm, arp) = pr.outer_wavespeeds(gamma_law_index);
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
        let p = u.to_primitive(gamma_law_index);
        println!("{:?} {:?}", primitive, p);
        assert!(f64::abs(primitive.velocity() - p.velocity()) < 1e-10);
        assert!(f64::abs(primitive.gas_pressure() - p.gas_pressure()) / primitive.gas_pressure() < 1e-10);
    }

    #[test]
    fn can_recover_primitive()
    {
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.0, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.2, 1.0));
        panic_unless_recovery_is_accurate(Primitive(1.0, 0.5, 1e-3));
        panic_unless_recovery_is_accurate(Primitive(1.0, 5.0, 1e+3));
    }
}
