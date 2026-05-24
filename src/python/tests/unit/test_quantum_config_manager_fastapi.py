from config import get_server_config
from quantum_calc.config_manager import QuantumCalculationConfigManager


def test_quantum_config_manager_has_no_framework_global_dependency():
    import quantum_calc.config_manager as module

    assert not hasattr(module, "current_app")


def test_quantum_config_manager_reads_server_config():
    manager = QuantumCalculationConfigManager()
    config = manager._get_config()

    server_config = get_server_config()
    expected = server_config.get("quantum_calculation_defaults", {})
    assert config["quantum_calculation_defaults"] == expected
    assert manager.get_memory_setting("DFT") == int(expected["memory_settings"]["DFT"])
    assert manager.get_config_info() == {
        "config_source": "ServerConfig",
        "config_available": True,
        "using_fallback": False,
    }
