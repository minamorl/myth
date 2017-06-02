extern crate num;
use std::iter::Sum;

#[derive(Debug)]
pub struct Matrix<T> {
    pub v: Vec<Vec<T>>,
    pub size: (usize, usize),
}

impl<T> Matrix<T> 
    where T: num::Signed + Clone + Copy + Sum, i32: Into<T> {
    pub fn new(v: Vec<Vec<T>>) -> Self {
        let size = (v.len(), v[0].len());
        Self {
            v: v,
            size: size,
        }
    }
    pub fn add(&self, vec: &Self) -> Result<Self, &str> {
        if self.size != vec.size {
            return Err("Matrix size is wrong");
        }
        Ok(Self::new((0..self.size.0).map(|i| {
            
            (0..self.size.1).map(|j| {
                self.v[i][j] + vec.v[i][j]
            }).collect()

        }).collect()))
    }
    pub fn transpose(&self) -> Self {
        Self::new((0..self.size.1).map(|j| {
            (0..self.size.0).map(|i| {
                self.v[i][j]
            }).collect()
        }).collect())
    }
    pub fn isSquare(&self) -> bool {
        self.size.0 == self.size.1
    }
    pub fn mul(&self, vec: &Self) -> Result<Self, &str> {
        if self.size.1 != vec.size.0 {
            Err("Wrong size")
        } else { 
            Ok(Self::new(
                (0..self.size.0).map(|i| {
                    (0..vec.size.1).map(|j| {
                        (0..vec.size.0).map(|k| self.v[i][k] * vec.v[k][j]).sum::<T>()
                    }).collect()
                }).collect()
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use matrix::Matrix;
    #[test]
    fn test_matrix_size() {
        let v = Matrix::new(
            vec![
                vec![1, 2, 3],
                vec![1, 2, 3],
            ]);
        assert_eq!(v.size, (2, 3));
    }
    #[test]
    fn test_matrix_add() {
        let v1 = Matrix::new(
            vec![
                vec![1, 2, 3],
                vec![1, 2, 3],
            ]);
        let v2 = Matrix::new(
            vec![
                vec![1, 2, 3],
                vec![1, 2, 3],
            ]);
        let v3 = v1.add(&v2).unwrap();
        assert_eq!(v3.v, vec![
            vec![2, 4, 6],
            vec![2, 4, 6],
        ]);
    }
    #[test]
    fn test_matrix_transpose() {
        let v1 = Matrix::new(
            vec![
                vec![1, 2, 3],
                vec![1, 2, 3],
            ]);
        let v3 = v1.transpose();
        assert_eq!(v3.v, vec![
            vec![1, 1],
            vec![2, 2],
            vec![3, 3]
        ]);
    }
    #[test]
    fn test_matrix_mul() {
        let v1 = Matrix::new(
            vec![
                vec![1, 2, 3],
                vec![4, 5, 6],
            ]);
        let v2 = v1.transpose();
        assert_eq!(v1.mul(&v2).unwrap().v, vec![
            vec![14, 32],
            vec![32, 77],
        ]);
    }
}