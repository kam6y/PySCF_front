-- PySCF_Front データベース初期スキーマ

-- インスタンス管理テーブル
CREATE TABLE instances (
    id VARCHAR(36) PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    status ENUM('draft', 'ready', 'running', 'completed', 'error') DEFAULT 'draft',
    user_id VARCHAR(36),
    project_id VARCHAR(36),
    INDEX idx_status (status),
    INDEX idx_user_id (user_id)
);

-- 分子情報テーブル
CREATE TABLE molecules (
    id VARCHAR(36) PRIMARY KEY,
    instance_id VARCHAR(36) NOT NULL,
    name VARCHAR(255),
    formula VARCHAR(100),
    molecular_weight DECIMAL(10,4),
    geometry_type ENUM('xyz', 'zmatrix') DEFAULT 'xyz',
    geometry_data JSON,
    charge INT DEFAULT 0,
    multiplicity INT DEFAULT 1,
    symmetry VARCHAR(10),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (instance_id) REFERENCES instances(id) ON DELETE CASCADE,
    INDEX idx_instance_id (instance_id)
);

-- 計算設定テーブル
CREATE TABLE calculations (
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
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (instance_id) REFERENCES instances(id) ON DELETE CASCADE,
    INDEX idx_instance_status (instance_id, status)
);

-- 計算結果テーブル
CREATE TABLE results (
    id VARCHAR(36) PRIMARY KEY,
    calculation_id VARCHAR(36) NOT NULL,
    result_type VARCHAR(50) NOT NULL,
    result_data JSON,
    file_path VARCHAR(500),
    file_size BIGINT,
    checksum VARCHAR(64),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (calculation_id) REFERENCES calculations(id) ON DELETE CASCADE,
    INDEX idx_calculation_type (calculation_id, result_type)
);

-- ジョブキュー管理テーブル
CREATE TABLE job_queue (
    id VARCHAR(36) PRIMARY KEY,
    calculation_id VARCHAR(36) NOT NULL,
    priority INT DEFAULT 5,
    status ENUM('waiting', 'running', 'completed', 'failed') DEFAULT 'waiting',
    assigned_worker VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP NULL,
    completed_at TIMESTAMP NULL,
    FOREIGN KEY (calculation_id) REFERENCES calculations(id) ON DELETE CASCADE,
    INDEX idx_status_priority (status, priority)
);

-- サンプルデータの挿入
INSERT INTO instances (id, name, description, status) VALUES 
('sample-instance-1', 'Sample Water Molecule', 'H2O molecule for testing', 'draft');

INSERT INTO molecules (id, instance_id, name, formula, molecular_weight, geometry_data, charge, multiplicity) VALUES 
('sample-mol-1', 'sample-instance-1', 'Water', 'H2O', 18.0153, 
 JSON_OBJECT(
   'atoms', JSON_ARRAY(
     JSON_OBJECT('symbol', 'O', 'x', 0.0, 'y', 0.0, 'z', 0.0),
     JSON_OBJECT('symbol', 'H', 'x', 0.7570, 'y', 0.5860, 'z', 0.0),
     JSON_OBJECT('symbol', 'H', 'x', -0.7570, 'y', 0.5860, 'z', 0.0)
   )
 ), 0, 1);