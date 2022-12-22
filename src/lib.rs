use itertools::Itertools;
use pyo3::{exceptions, prelude::*};
use seal::pair::{AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, Step};
use std::collections::HashMap;

#[pyclass]
struct AlignmentResult {
    #[pyo3(get)]
    steps: Vec<u8>,
    #[pyo3(get)]
    score: isize,
}

fn matrix_match_fn<'a>(
    similarity_matrix: Option<HashMap<(&'a str, &'a str), isize>>,
) -> impl Fn(&'a str, &'a str) -> isize {
    move |x: &'a str, y: &'a str| -> isize {
        if let Some(matrix) = &similarity_matrix {
            match matrix.get(&(x, y)) {
                Some(score) => *score,
                None => match matrix.get(&(y, x)) {
                    Some(score) => *score,

                    None => {
                        if x == y {
                            1
                        } else {
                            -1
                        }
                    }
                },
            }
        } else {
            if x == y {
                1
            } else {
                -1
            }
        }
    }
}

/// Finds alignment similarity between two sequences
#[pyfunction]
fn align(
    _py: Python,
    a: Vec<&str>,
    b: Vec<&str>,
    similarity_matrix: Option<HashMap<(&str, &str), isize>>,
) -> PyResult<AlignmentResult> {
    let needleman_wunsch = NeedlemanWunsch::new(-1, -1, -1);

    let match_fn = matrix_match_fn(similarity_matrix);
    let alignment_set: Result<AlignmentSet<InMemoryAlignmentMatrix>, _> =
        AlignmentSet::new(a.len(), b.len(), needleman_wunsch, |x, y| {
            match_fn(a[x], b[y])
        });

    match alignment_set {
        Ok(ref alignment_set) => {
            let global_alignment = alignment_set.global_alignment();
            // trace(&a, &b, &global_alignment);
            // Ok(global_alignment.score())
            Ok(AlignmentResult {
                steps: global_alignment
                    .steps()
                    .map(|step| match step {
                        Step::Align { x: _, y: _ } => 0,
                        Step::Delete { x: _ } => 1,
                        Step::Insert { y: _ } => 2,
                    })
                    .collect(),
                score: global_alignment.score(),
            })
        }
        Err(error) => Err(exceptions::PyValueError::new_err(error)),
    }
}

/// Takes a set of sequences and finds all pair-wise similarity distances
#[pyfunction]
fn align_set(
    _py: Python,
    seqs: Vec<Vec<&str>>,
    similarity_matrix: Option<HashMap<(&str, &str), isize>>,
) -> PyResult<HashMap<(usize, usize), AlignmentResult>> {
    let needleman_wunsch = NeedlemanWunsch::new(-1, -1, -1);
    let match_fn = matrix_match_fn(similarity_matrix);
    let mut scores = HashMap::new();
    let len = seqs.len();

    for pair in (0..len).combinations(2) {
        let x = pair[0];
        let y = pair[1];
        let a = &seqs[x];
        let b = &seqs[y];

        let alignment_set: Result<AlignmentSet<InMemoryAlignmentMatrix>, _> =
            AlignmentSet::new(a.len(), b.len(), needleman_wunsch.clone(), |x, y| {
                match_fn(a[x], b[y])
            });

        match alignment_set {
            Ok(alignment_set) => {
                let global_alignment = alignment_set.global_alignment();
                scores.insert(
                    (x, y),
                    // global_alignment.score(),
                    AlignmentResult {
                        steps: global_alignment
                            .steps()
                            .map(|step| match step {
                                Step::Align { x: _, y: _ } => 0,
                                Step::Delete { x: _ } => 1,
                                Step::Insert { y: _ } => 2,
                            })
                            .collect(),
                        score: global_alignment.score(),
                    },
                )
            }
            Err(error) => return Err(exceptions::PyValueError::new_err(error)),
        };
    }

    Ok(scores)
}

/// A Python module implemented in Rust.
#[pymodule]
fn sequences(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align, m)?)?;
    m.add_function(wrap_pyfunction!(align_set, m)?)?;
    m.add_class::<AlignmentResult>()?;
    Ok(())
}
