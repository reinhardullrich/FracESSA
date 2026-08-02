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
    candidate_count INTEGER CHECK (candidate_count >= 0),
    ess_count INTEGER CHECK (ess_count BETWEEN 0 AND candidate_count),
    candidate_structure TEXT CHECK (
        candidate_structure IS NULL OR (
            json_valid(candidate_structure) AND
            json_type(candidate_structure) = 'object'
        )
    ),
    ess_structure TEXT CHECK (
        ess_structure IS NULL OR (
            json_valid(ess_structure) AND
            json_type(ess_structure) = 'object'
        )
    ),
    origin TEXT NOT NULL CHECK (length(origin) > 0),
    tags TEXT NOT NULL DEFAULT '[]' CHECK (
        json_valid(tags) AND
        json_type(tags) = 'array'
    ),
    name TEXT CHECK (name IS NULL OR length(name) > 0),
    family TEXT CHECK (family IS NULL OR length(family) > 0),
    subfamily TEXT CHECK (subfamily IS NULL OR length(subfamily) > 0),
    description TEXT CHECK (description IS NULL OR length(description) > 0),
    source_url TEXT CHECK (source_url IS NULL OR length(source_url) > 0),
    original_format TEXT CHECK (
        original_format IS NULL OR length(original_format) > 0
    ),
    original_id TEXT CHECK (original_id IS NULL OR length(original_id) > 0),
    created_at TEXT CHECK (created_at IS NULL OR length(created_at) > 0),
    CHECK (
        (candidate_count IS NULL AND ess_count IS NULL AND
         candidate_structure IS NULL AND ess_structure IS NULL) OR
        (candidate_count IS NOT NULL AND ess_count IS NOT NULL AND
         candidate_structure IS NOT NULL AND ess_structure IS NOT NULL)
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
    multiplier INTEGER CHECK (multiplier IS NULL OR multiplier >= 1),
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

CREATE TABLE timings (
    session TEXT NOT NULL CHECK (length(session) > 0),
    recorded_at TEXT NOT NULL CHECK (length(recorded_at) > 0),
    machine TEXT NOT NULL CHECK (length(machine) > 0),
    cpu_id INTEGER NOT NULL CHECK (cpu_id >= 0),
    comment TEXT NOT NULL DEFAULT '',
    build_label TEXT NOT NULL CHECK (length(build_label) > 0),
    source_ref TEXT NOT NULL CHECK (length(source_ref) > 0),
    revision TEXT NOT NULL CHECK (length(revision) > 0),
    binary_sha256 TEXT NOT NULL CHECK (length(binary_sha256) = 64),
    backend TEXT NOT NULL CHECK (backend IN ('pybind', 'cli')),
    mode TEXT NOT NULL CHECK (mode IN ('fast', 'safe', 'unsafe', 'exact')),
    matrix_id INTEGER NOT NULL,
    target_ns INTEGER NOT NULL CHECK (target_ns > 0),
    iterations INTEGER NOT NULL CHECK (iterations > 0),
    measured_wall_ns INTEGER NOT NULL CHECK (measured_wall_ns > 0),
    elapsed_ns INTEGER NOT NULL CHECK (elapsed_ns >= 0),
    ess_count INTEGER NOT NULL CHECK (ess_count >= 0),
    PRIMARY KEY (session, build_label, mode, matrix_id),
    FOREIGN KEY (matrix_id) REFERENCES matrices(matrix_id) ON DELETE CASCADE
) STRICT;

CREATE VIEW timing_overview AS
SELECT
    m.matrix_id,
    m.dimension,
    m.is_cs,
    m.name,
    t.elapsed_ns AS calibration_ns,
    t.elapsed_ns / 1000000000.0 AS calibration_seconds,
    t.build_label,
    t.mode,
    t.session,
    t.recorded_at,
    t.target_ns,
    t.iterations,
    t.measured_wall_ns,
    t.ess_count AS measured_ess_count,
    m.ess_count AS expected_ess_count,
    t.machine,
    t.cpu_id,
    t.backend,
    t.source_ref,
    t.revision,
    t.binary_sha256,
    t.comment,
    m.size_class,
    m.candidate_count,
    m.candidate_structure,
    m.ess_structure,
    m.family,
    m.subfamily,
    m.origin,
    m.tags,
    m.description,
    m.source_url,
    m.original_format,
    m.original_id,
    m.created_at,
    m.matrix
FROM matrices AS m
LEFT JOIN timings AS t USING (matrix_id)
ORDER BY m.dimension, m.is_cs, m.matrix_id,
         t.recorded_at, t.build_label, t.mode;
