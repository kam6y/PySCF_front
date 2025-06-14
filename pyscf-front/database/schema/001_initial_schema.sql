-- PySCF_Front Database Schema v1.0
-- Initial schema for quantum chemistry calculations

-- Instances table - Main calculation instances
CREATE TABLE IF NOT EXISTS instances (
    id VARCHAR(36) PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    status ENUM('draft', 'ready', 'running', 'completed', 'error', 'cancelled') DEFAULT 'draft',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    user_id VARCHAR(36),
    project_id VARCHAR(36),
    INDEX idx_status (status),
    INDEX idx_user_id (user_id),
    INDEX idx_created_at (created_at)
);

-- Molecules table - Molecular structure data
CREATE TABLE IF NOT EXISTS molecules (
    id VARCHAR(36) PRIMARY KEY,
    instance_id VARCHAR(36) NOT NULL,
    name VARCHAR(255),
    formula VARCHAR(100),
    molecular_weight DECIMAL(10,4),
    geometry_type ENUM('xyz', 'zmatrix') DEFAULT 'xyz',
    geometry_data JSON NOT NULL,
    charge INT DEFAULT 0,
    multiplicity INT DEFAULT 1,
    symmetry VARCHAR(10),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (instance_id) REFERENCES instances(id) ON DELETE CASCADE,
    INDEX idx_instance_id (instance_id),
    INDEX idx_formula (formula)
);

-- Calculations table - Calculation parameters and status
CREATE TABLE IF NOT EXISTS calculations (
    id VARCHAR(36) PRIMARY KEY,
    instance_id VARCHAR(36) NOT NULL,
    method VARCHAR(50) NOT NULL,
    basis_set VARCHAR(50) NOT NULL,
    parameters JSON,
    convergence_criteria JSON,
    max_iterations INT DEFAULT 100,
    start_time TIMESTAMP NULL,
    end_time TIMESTAMP NULL,
    status ENUM('pending', 'running', 'completed', 'failed', 'cancelled') DEFAULT 'pending',
    error_message TEXT,
    calculation_time_seconds DECIMAL(10,3),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (instance_id) REFERENCES instances(id) ON DELETE CASCADE,
    INDEX idx_instance_status (instance_id, status),
    INDEX idx_method_basis (method, basis_set),
    INDEX idx_status (status)
);

-- Results table - Calculation results and output data
CREATE TABLE IF NOT EXISTS results (
    id VARCHAR(36) PRIMARY KEY,
    calculation_id VARCHAR(36) NOT NULL,
    result_type VARCHAR(50) NOT NULL, -- 'energy', 'orbitals', 'properties', 'geometry', etc.
    result_data JSON NOT NULL,
    file_path VARCHAR(500),
    file_size BIGINT,
    checksum VARCHAR(64),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (calculation_id) REFERENCES calculations(id) ON DELETE CASCADE,
    INDEX idx_calculation_type (calculation_id, result_type),
    INDEX idx_result_type (result_type)
);

-- Job queue table - Calculation job management
CREATE TABLE IF NOT EXISTS job_queue (
    id VARCHAR(36) PRIMARY KEY,
    calculation_id VARCHAR(36) NOT NULL,
    priority INT DEFAULT 5,
    status ENUM('waiting', 'running', 'completed', 'failed') DEFAULT 'waiting',
    assigned_worker VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP NULL,
    completed_at TIMESTAMP NULL,
    error_message TEXT,
    retry_count INT DEFAULT 0,
    max_retries INT DEFAULT 3,
    FOREIGN KEY (calculation_id) REFERENCES calculations(id) ON DELETE CASCADE,
    INDEX idx_status_priority (status, priority),
    INDEX idx_assigned_worker (assigned_worker),
    INDEX idx_created_at (created_at)
);

-- Projects table (optional) - For organizing calculations
CREATE TABLE IF NOT EXISTS projects (
    id VARCHAR(36) PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    created_by VARCHAR(36),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    INDEX idx_created_by (created_by),
    INDEX idx_name (name)
);

-- System configuration table
CREATE TABLE IF NOT EXISTS system_config (
    id INT AUTO_INCREMENT PRIMARY KEY,
    config_key VARCHAR(100) UNIQUE NOT NULL,
    config_value TEXT,
    description TEXT,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    INDEX idx_config_key (config_key)
);

-- Insert default system configuration
INSERT INTO system_config (config_key, config_value, description) VALUES
('max_concurrent_jobs', '4', 'Maximum number of concurrent calculation jobs'),
('default_timeout_minutes', '60', 'Default calculation timeout in minutes'),
('cleanup_completed_jobs_days', '30', 'Days to keep completed job records'),
('enable_gpu_acceleration', 'false', 'Enable GPU acceleration if available'),
('log_level', 'INFO', 'System logging level')
ON DUPLICATE KEY UPDATE config_value = VALUES(config_value);

-- Sample data for development
INSERT INTO instances (id, name, description, status) VALUES
('sample-instance-1', 'Water Molecule Test', 'Basic water molecule calculation for testing', 'draft'),
('sample-instance-2', 'Methane Molecule Test', 'Methane molecule optimization test', 'draft')
ON DUPLICATE KEY UPDATE name = VALUES(name);

INSERT INTO molecules (id, instance_id, name, formula, molecular_weight, geometry_data, charge, multiplicity) VALUES
(
    'molecule-water-1', 
    'sample-instance-1', 
    'Water', 
    'H2O', 
    18.0153,
    JSON_OBJECT(
        'atoms', JSON_ARRAY(
            JSON_OBJECT('symbol', 'O', 'x', 0.000000, 'y', 0.000000, 'z', 0.000000),
            JSON_OBJECT('symbol', 'H', 'x', 0.757000, 'y', 0.586000, 'z', 0.000000),
            JSON_OBJECT('symbol', 'H', 'x', -0.757000, 'y', 0.586000, 'z', 0.000000)
        ),
        'atom_count', 3
    ),
    0,
    1
),
(
    'molecule-methane-1', 
    'sample-instance-2', 
    'Methane', 
    'CH4', 
    16.0425,
    JSON_OBJECT(
        'atoms', JSON_ARRAY(
            JSON_OBJECT('symbol', 'C', 'x', 0.000000, 'y', 0.000000, 'z', 0.000000),
            JSON_OBJECT('symbol', 'H', 'x', 1.089000, 'y', 0.000000, 'z', 0.000000),
            JSON_OBJECT('symbol', 'H', 'x', -0.363000, 'y', 1.026000, 'z', 0.000000),
            JSON_OBJECT('symbol', 'H', 'x', -0.363000, 'y', -0.513000, 'z', 0.889000),
            JSON_OBJECT('symbol', 'H', 'x', -0.363000, 'y', -0.513000, 'z', -0.889000)
        ),
        'atom_count', 5
    ),
    0,
    1
)
ON DUPLICATE KEY UPDATE name = VALUES(name);