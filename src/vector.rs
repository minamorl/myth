extern crate num;
use std::ops::*;

pub struct Vector<T> {
    pub v: Vec<T>,
    pub dim: usize,
}
impl<T> Vector<T> 
    where T: num::Signed + Clone + Copy, i32: Into<T> {
    pub fn new(v: Vec<T>) -> Self {
        let dim = v.len();
        Self {
            v: v,
            dim: dim,
        }
    }
    pub fn add(&self, vec: &Self) -> Result<Self, &str> {
        if self.dim != vec.dim {
            Err("Vec size is wrong")
        } else {
            Ok(Vector::new((0..self.dim).map(|i| self.v[i] + vec.v[i]).collect()))
        }
    }
    pub fn inverse(&self) -> Self {
        let minus: T = -1i32.into();
        self.scalar(minus)
    }
    pub fn dot(&self, vec: &Self) -> Result<T, &str> {
        if self.dim != vec.dim {
            Err("Vec size is wrong")
        } else {
            Ok((0..self.dim).map(|i| self.v[i] * vec.v[i]).fold(0.into(), |sum, val| sum + val))
        }
    }
    pub fn cross(&self, vec: &Self) -> Result<Self, &str> {
        let len1 = self.dim;
        let len2 = vec.dim;
        if len1 != len2 || len1 < 2 || len1 > 3 || len2 < 2 || len2 > 3 {
            return Err("Vec size should be 2 or 3");
        }
        if len1 == 2 {
            Ok(Self::new(vec![self.v[0] * vec.v[1] - self.v[1] * vec.v[0]])) 
        } else {
            Ok(Self::new(vec![
                self.v[1] * vec.v[2] - self.v[2] * vec.v[1],
                self.v[2] * vec.v[0] - self.v[0] * vec.v[2],   
                self.v[0] * vec.v[1] - self.v[1] * vec.v[0],
            ]))
        }
    }
    pub fn scalar<U>(&self, n: U) -> Self where U: num::Signed + Into<T> + Copy + Clone { 
        Self::new(
            self.v.clone().iter().map(|x| n.into() * *(x)).collect()
        )
    }
}

#[cfg(test)]
mod tests {
    use vector::Vector;

    #[test]
    fn test_vector_add() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let v2 = Vector::new(vec![1, 2, 3]);
        let v3 = v1.add(&v2).unwrap();
        assert_eq!(v3.v, vec![2, 4, 6]);
    }

    #[test]
    fn test_vector_dot() {
        let v1 = Vector::new(vec![1.0, 2.0, 3.0]);
        let v2 = Vector::new(vec![1.0, 2.0, 3.0]);
        let r = v1.dot(&v2).unwrap();
        assert_eq!(r, 14.0);
    }

    #[test]
    fn test_vector_inverse() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let v2 = Vector::new(vec![-1, -2, -3]);
        let r = v1.inverse();
        assert_eq!(r.v, v2.v);
    }
    #[test]
    fn test_vector_cross() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let v2 = Vector::new(vec![1, 2, 3]);
        let r = v1.cross(&v2).unwrap();
        assert_eq!(r.v, vec![2 * 3 - 3 * 2, 3 * 1 - 1 * 3, 1 * 2 - 2 * 1]);
    }

    #[test]
    fn test_vector_scalar() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let r = v1.scalar(3);
        assert_eq!(r.v, vec![3, 6, 9]);
    }
}