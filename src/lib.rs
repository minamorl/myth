extern crate num;

use std::ops::*;

struct Vector<T> {
    v: Vec<T>,
    dim: usize,
}
impl<T> Vector<T> 
    where T: num::Signed + Clone + Copy, i32: Into<T> {
    pub fn new(v: Vec<T>) -> Self {
        let dim = v.len();
        Vector {
            v: v,
            dim: dim,
        }
    }
    fn add(&self, vec: &Self) -> Self {
        if self.dim != vec.dim {
            panic!("Vec size is wrong");
        }
        let mut r = vec![];
        for i in 0..self.dim {
            r.push(self.v[i] + vec.v[i]);
        }
        Vector::new(r)
    }
    fn inverse(&self) -> Self {
        let minus: T = -1i32.into();
        self.scalar(minus)
    }
    fn dot(&self, vec: &Self) -> T {
        if self.dim != vec.dim {
            panic!("Vec size is wrong");
        }
        (0..self.dim).map(|i| self.v[i] * vec.v[i]).fold(0.into(), |sum, val| sum + val)
    }
    fn cross(&self, vec: &Self) -> Self {
        let len1 = self.dim;
        let len2 = vec.dim;
        if len1 != len2 || len1 < 2 || len1 > 3 || len2 < 2 || len2 > 3 {
            panic!("Vec size should be 2 or 3");
        }
        if len1 == 2 {
            Self::new(vec![self.v[0] * vec.v[1] - self.v[1] * vec.v[0]]) 
        } else {
            Self::new(vec![
                self.v[1] * vec.v[2] - self.v[2] * vec.v[1],
                self.v[2] * vec.v[0] - self.v[0] * vec.v[2],   
                self.v[0] * vec.v[1] - self.v[1] * vec.v[0],
            ])
        }
    }
    fn scalar<U>(&self, n: U) -> Self where U: num::Signed + Into<T> + Copy + Clone { 
        Self::new(
          self.v.clone().iter().map(|x| n.into() * *(x)).collect()
        )
    }
}

#[cfg(test)]
mod tests {
    use Vector;
    #[test]
    fn test_vector_add() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let v2 = Vector::new(vec![1, 2, 3]);
        let v3 = v1.add(&v2);
        assert_eq!(v3.v, vec![2, 4, 6]);
    }

    #[test]
    fn test_vector_dot() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let v2 = Vector::new(vec![1, 2, 3]);
        let r = v1.dot(&v2);
        assert_eq!(r, 14);
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
        let r = v1.cross(&v2);
        assert_eq!(r.v, vec![2 * 3 - 3 * 2, 3 * 1 - 1 * 3, 1 * 2 - 2 * 1]);
    }

    #[test]
    fn test_vector_scalar() {
        let v1 = Vector::new(vec![1, 2, 3]);
        let r = v1.scalar(3);
        assert_eq!(r.v, vec![3, 6, 9]);
    }
}