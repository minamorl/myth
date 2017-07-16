extern crate num;
use std::iter::{Sum};
use std::ops::Mul;
use std::cmp::{Ord};
#[derive(Debug)]
pub struct Matrix<T> {
    pub v: Vec<Vec<T>>,
    pub size: (usize, usize),
}

#[derive(Debug)]
pub struct LUPDecomposed<T> {
    pub P: Matrix<T>,
    pub L: Matrix<T>,
    pub U: Matrix<T>,
}

impl<T> Matrix<T> 
    where T: num::Signed + Clone + Copy + Sum + Mul + PartialOrd, i32: Into<T> {
    pub fn new(v: Vec<Vec<T>>) -> Self {
        let size = (v.len(), v[0].len());
        Self {
            v: v,
            size: size,
        }
    }
    pub fn zeroes(n: usize, m: usize) -> Self {
        Self::new(
            vec![vec![0.into(); m]; n]
        )
    }
    pub fn diag(v: Vec<T>) -> Self {
        Self::new(
            (0..v.len()).map(|i| {
                (0..v.len()).map(|j| {
                    if i == j {
                        v[i]
                    } else {
                        0.into()
                    }
                }).collect()
            }).collect()
        )
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
    pub fn scalar<U: Into<T> + Copy + Clone>(&self, n: U) -> Self {
        Self::new(
            self.v.clone().iter().map(|l| {
                l.iter().map(|x| (*x) * n.into()).collect()
            }).collect()
        )
    }
    pub fn pivot(&self) -> Self {
        let mut pivot = Self::diag(vec![1.into(); self.size.0]);

        for j in 0..self.size.0 {
            for i in j..self.size.0 {
                if self.v[i][j].abs() > self.v[j][j].abs() {
                    pivot.v.swap(i, j);
                }
           }
        }
        pivot

    }
    /// Provides a result from LU decomposition
    pub fn decompose(&self) -> LUPDecomposed<T> {
        let mut L = Matrix::diag(vec![1.into(); self.size.0]);
        let mut U = Matrix::zeroes(self.size.0, self.size.0);
        let P = self.pivot();
        let PA = P.mul(&self).unwrap();
        for j in 0..self.v.len() {
            for i in 0..j + 1 {
                U.v[i][j] = PA.v[i][j] - (0..i).map(|k| U.v[k][j] * L.v[i][k]).sum()
            }
            for i in j..self.size.0 {
                L.v[i][j] = (PA.v[i][j] - (0..j).map(|k| L.v[i][k] * U.v[k][j]).sum()) / U.v[j][j]
            }
        }
        LUPDecomposed {
            P: P,
            L: L,
            U: U,
        }
    }
}

impl<T> LUPDecomposed<T> 
    where T: num::Signed + Clone + Copy + Sum + Mul + PartialOrd, i32: Into<T>{

    /// Determinant from PLU
    pub fn determinant(&self) -> T {
        let mut flag: T = -1.into();
        let mut lsum: T = 1.into();
        let mut usum: T = 1.into();
        for i in 0..self.P.size.0 {
            flag = flag * (if self.P.v[i][i] == 0.into() { 1.into() } else { -1.into() });
            lsum = lsum * self.L.v[i][i];
            usum = usum * self.U.v[i][i];
        }
        flag * lsum * usum
    }
}

#[cfg(test)]
mod tests {
    use matrix::Matrix;
    use matrix::LUPDecomposed;

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
    #[test]
    fn test_matrix_diag() {
        let v1 = 
            vec![
                vec![1, 0, 0],
                vec![0, 1, 0],
                vec![0, 0, 1],
            ];
        assert_eq!(Matrix::diag(vec![1; 3]).v, v1);
    }
    #[test]
    fn test_matrix_scalar() {
        assert_eq!(Matrix::diag(vec![1; 3]).scalar(2).v, Matrix::diag(vec![2; 3]).v);
    }
    #[test]
    fn test_matrix_decompose() {
        let m = Matrix::new(
            vec![
                vec![7.0, 3.0, -1.0, 2.0],
                vec![3.0, 8.0, 1.0, -4.0],
                vec![-1.0, 1.0, 4.0, -1.0],
                vec![2.0, -4.0, -1.0, 6.0],
            ]
        );
        let decomposed = m.decompose();

        assert_eq!(decomposed.L.v, vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.42857142857142855, 1.0, 0.0, 0.0],
            vec![-0.14285714285714285, 0.2127659574468085, 1.0, 0.0],
            vec![0.2857142857142857, -0.7234042553191489, 0.0898203592814371, 1.0]
        ]);

        assert_eq!(decomposed.U.v, vec![
            vec![7.0, 3.0, -1.0, 2.0], 
            vec![0.0, 6.714285714285714, 1.4285714285714286, -4.857142857142857],
            vec![0.0, 0.0, 3.5531914893617023, 0.31914893617021267],
            vec![0.0, 0.0, 0.0, 1.88622754491018],
        ]);
    }
    #[test]
    fn test_LUPDecomposed_determinant() {
        let m = Matrix::new(
            vec![
                vec![3.0, 1.0, 1.0, 2.0],
                vec![5.0, 1.0, 3.0, 4.0],
                vec![2.0, 0.0, 1.0, 0.0],
                vec![1.0, 3.0, 2.0, 1.0],
            ]
        );
        assert_eq!(m.decompose().determinant().round(), -22.0)
    }
}