use pyo3::{exceptions, prelude::*};
use seal::pair::{Alignment, AlignmentSet, InMemoryAlignmentMatrix, NeedlemanWunsch, Step};
use std::collections::HashMap;

type SimilarityMatrix<'a> = HashMap<(&'a str, &'a str), isize>;

fn trace<'a, T: ToString + Copy>(
    x_seq: &'a Vec<T>,
    y_seq: &'a Vec<T>,
    alignment: &'a Alignment,
) -> impl Iterator<Item = (String, String)> + 'a {
    alignment.steps().map(move |step| match step {
        Step::Align { x, y } => (x_seq[x].to_string(), y_seq[y].to_string()),
        Step::Delete { x } => (x_seq[x].to_string(), String::from("-")),
        Step::Insert { y } => (String::from("-"), y_seq[y].to_string()),
    })
}

#[pyclass]
struct AlignmentResult {
    #[pyo3(get)]
    alignments: Vec<(String, String)>,
    #[pyo3(get)]
    alignment_score: isize,
    #[pyo3(get)]
    similarity_score: f64,
}

struct Scorer<'a> {
    matrix: &'a SimilarityMatrix<'a>,
    match_score: isize,
    mismatch_score: isize,
}

impl Scorer<'_> {
    fn compare(&self, x: &str, y: &str) -> isize {
        match self.matrix.get(&(x, y)) {
            Some(score) => *score,
            None => match self.matrix.get(&(y, x)) {
                Some(score) => *score,

                None => {
                    if x == y {
                        self.match_score
                    } else {
                        self.mismatch_score
                    }
                }
            },
        }
    }

    fn similarity_score(&self, x_seq: &Vec<&str>, y_seq: &Vec<&str>, alignment: &Alignment) -> f64 {
        let (dis_correct, num_correct): (i32, u32) =
            alignment.steps().fold((0, 0), |(dc, nc), step| match step {
                Step::Align { x, y } => {
                    if x_seq[x] == y_seq[y] {
                        (dc + self.compare(&x_seq[x], &y_seq[y]) as i32, nc + 1)
                    } else {
                        (dc, nc)
                    }
                }
                _ => (dc, nc),
            });

        if num_correct == 0 {
            return -1f64;
        }

        let dis = alignment.score() as i32;

        let sim_align = match dis_correct {
            0 => 0f64,
            _ => f64::from(dis) / f64::from(dis_correct),
        };

        let sim_significance = f64::from(num_correct) / f64::from(alignment.len() as i32);

        sim_align * sim_significance
    }
}

/// Finds alignment similarity between two sequences
#[pyfunction(match_score=1, mismatch_score=-1, gap_score=-1)]
fn align(
    _py: Python,
    a: Vec<&str>,
    b: Vec<&str>,
    match_score: isize,
    mismatch_score: isize,
    gap_score: isize,
    similarity_matrix: Option<SimilarityMatrix>,
) -> PyResult<AlignmentResult> {
    let needleman_wunsch = NeedlemanWunsch::new(mismatch_score, gap_score, gap_score);
    let scorer = Scorer {
        matrix: &similarity_matrix.unwrap_or(HashMap::new()),
        match_score,
        mismatch_score,
    };

    let alignment_set: Result<AlignmentSet<InMemoryAlignmentMatrix>, _> =
        AlignmentSet::new(a.len(), b.len(), needleman_wunsch, |x, y| {
            scorer.compare(a[x], b[y])
        });

    match alignment_set {
        Ok(ref alignment_set) => {
            let global_alignment = alignment_set.global_alignment();
            Ok(AlignmentResult {
                alignments: trace(&a, &b, &global_alignment).collect(),
                alignment_score: global_alignment.score(),
                similarity_score: scorer.similarity_score(&a, &b, &global_alignment),
            })
        }
        Err(error) => Err(exceptions::PyValueError::new_err(error)),
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn sequences(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align, m)?)?;
    m.add_class::<AlignmentResult>()?;
    Ok(())
}
