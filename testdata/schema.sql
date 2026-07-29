PRAGMA foreign_keys = ON;

CREATE TABLE matrices (
    matrix_id INTEGER PRIMARY KEY,
    dimension INTEGER NOT NULL CHECK (dimension BETWEEN 1 AND 63),
    size_class TEXT NOT NULL CHECK (
        (dimension BETWEEN 1 AND 8 AND size_class = 'small') OR
        (dimension BETWEEN 9 AND 16 AND size_class = 'medium') OR
        (dimension >= 17 AND size_class = 'large')
    ),
    is_cs INTEGER NOT NULL CHECK (is_cs IN (0, 1)),
    matrix TEXT NOT NULL CHECK (length(matrix) > 0),
    candidate_count INTEGER NOT NULL CHECK (candidate_count >= 0),
    ess_count INTEGER NOT NULL CHECK (ess_count BETWEEN 0 AND candidate_count),
    candidate_structure TEXT NOT NULL CHECK (
        json_valid(candidate_structure) AND
        json_type(candidate_structure) = 'object'
    ),
    ess_structure TEXT NOT NULL CHECK (
        json_valid(ess_structure) AND
        json_type(ess_structure) = 'object'
    ),
    origin TEXT NOT NULL CHECK (length(origin) > 0),
    tags TEXT NOT NULL DEFAULT '[]' CHECK (
        json_valid(tags) AND
        json_type(tags) = 'array'
    )
) STRICT;

CREATE TABLE candidates (
    matrix_id INTEGER NOT NULL,
    candidate_id INTEGER NOT NULL CHECK (candidate_id > 0),
    vector TEXT NOT NULL,
    support INTEGER NOT NULL CHECK (support > 0),
    support_size INTEGER NOT NULL CHECK (support_size BETWEEN 1 AND 63),
    extended_support INTEGER NOT NULL CHECK (extended_support > 0),
    extended_support_size INTEGER NOT NULL CHECK (
        extended_support_size BETWEEN support_size AND 63
    ),
    shift_reference INTEGER NOT NULL CHECK (shift_reference >= 0),
    is_ess INTEGER NOT NULL CHECK (is_ess IN (0, 1)),
    reason_ess TEXT NOT NULL,
    payoff TEXT NOT NULL,
    payoff_double TEXT NOT NULL,
    PRIMARY KEY (matrix_id, candidate_id),
    UNIQUE (matrix_id, support),
    FOREIGN KEY (matrix_id) REFERENCES matrices(matrix_id) ON DELETE CASCADE
) STRICT;

CREATE INDEX candidates_by_matrix_and_ess
    ON candidates(matrix_id, is_ess, support_size);
