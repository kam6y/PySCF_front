"""System resource manager for monitoring CPU and memory usage."""

import os
import logging
import multiprocessing
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime, timezone

try:
    import psutil
except ImportError:
    psutil = None
    
logger = logging.getLogger(__name__)


@dataclass
class SystemResourceInfo:
    """System resource information."""
    total_cpu_cores: int
    total_memory_mb: int
    available_memory_mb: int
    cpu_usage_percent: float
    memory_usage_percent: float
    timestamp: datetime


@dataclass
class CalculationResourceUsage:
    """Resource usage for a specific calculation."""
    calculation_id: str
    cpu_cores: int
    memory_mb: int
    estimated_memory_mb: int  # Estimated memory usage based on calculation type
    timestamp: datetime


@dataclass
class ResourceConstraints:
    """Resource usage constraints for the system."""
    max_cpu_utilization_percent: float = 80.0
    max_memory_utilization_percent: float = 80.0
    system_total_cores: int = 0
    system_total_memory_mb: int = 0


class SystemResourceManager:
    """Manager for system resource monitoring and allocation."""
    
    def __init__(self):
        """Initialize the system resource manager."""
        self._active_calculations: Dict[str, CalculationResourceUsage] = {}
        self._resource_constraints = ResourceConstraints()
        
        if psutil is None:
            logger.warning("psutil is not available. Resource monitoring will be limited.")
            # Fallback to basic multiprocessing info
            self._resource_constraints.system_total_cores = multiprocessing.cpu_count()
            # Estimate 4GB as fallback (conservative)
            self._resource_constraints.system_total_memory_mb = 4096
        else:
            # Get system information
            self._resource_constraints.system_total_cores = psutil.cpu_count(logical=True)
            self._resource_constraints.system_total_memory_mb = int(psutil.virtual_memory().total / (1024 * 1024))
            
        logger.info(f"Initialized SystemResourceManager - CPU cores: {self._resource_constraints.system_total_cores}, "
                   f"Total memory: {self._resource_constraints.system_total_memory_mb} MB")
    
    def get_system_info(self) -> SystemResourceInfo:
        """Get current system resource information."""
        if psutil is None:
            # Return conservative estimates when psutil is not available
            return SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=self._resource_constraints.system_total_memory_mb,
                available_memory_mb=int(self._resource_constraints.system_total_memory_mb * 0.5),  # Assume 50% available
                cpu_usage_percent=50.0,  # Conservative estimate
                memory_usage_percent=50.0,  # Conservative estimate
                timestamp=datetime.now(timezone.utc)
            )
        
        try:
            memory = psutil.virtual_memory()
            cpu_percent = psutil.cpu_percent(interval=0.1)
            
            return SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=int(memory.total / (1024 * 1024)),
                available_memory_mb=int(memory.available / (1024 * 1024)),
                cpu_usage_percent=cpu_percent,
                memory_usage_percent=memory.percent,
                timestamp=datetime.now(timezone.utc)
            )
        except Exception as e:
            logger.error(f"Error getting system information: {e}")
            # Return fallback values
            return SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=self._resource_constraints.system_total_memory_mb,
                available_memory_mb=int(self._resource_constraints.system_total_memory_mb * 0.5),
                cpu_usage_percent=0.0,
                memory_usage_percent=0.0,
                timestamp=datetime.now(timezone.utc)
            )
    
    def get_resource_constraints(self) -> ResourceConstraints:
        """Get current resource constraints."""
        return self._resource_constraints
    
    def update_resource_constraints(self, 
                                  max_cpu_utilization_percent: Optional[float] = None,
                                  max_memory_utilization_percent: Optional[float] = None) -> None:
        """Update resource constraints."""
        if max_cpu_utilization_percent is not None:
            self._resource_constraints.max_cpu_utilization_percent = max(10.0, min(100.0, max_cpu_utilization_percent))
        
        if max_memory_utilization_percent is not None:
            self._resource_constraints.max_memory_utilization_percent = max(10.0, min(100.0, max_memory_utilization_percent))
            
        logger.info(f"Updated resource constraints - CPU: {self._resource_constraints.max_cpu_utilization_percent}%, "
                   f"Memory: {self._resource_constraints.max_memory_utilization_percent}%")
    
    def register_calculation(self, calculation_id: str, cpu_cores: int, memory_mb: int,
                           calculation_method: str = "DFT") -> None:
        """Register a calculation with its resource requirements."""
        # Estimate memory usage based on calculation method
        estimated_memory = self._estimate_memory_usage(memory_mb, calculation_method, cpu_cores)
        
        calculation_usage = CalculationResourceUsage(
            calculation_id=calculation_id,
            cpu_cores=cpu_cores,
            memory_mb=memory_mb,
            estimated_memory_mb=estimated_memory,
            timestamp=datetime.now(timezone.utc)
        )
        
        self._active_calculations[calculation_id] = calculation_usage
        logger.info(f"Registered calculation {calculation_id} - CPU cores: {cpu_cores}, "
                   f"Memory: {memory_mb} MB (estimated: {estimated_memory} MB)")
    
    def unregister_calculation(self, calculation_id: str) -> None:
        """Unregister a calculation when it completes."""
        if calculation_id in self._active_calculations:
            calculation = self._active_calculations.pop(calculation_id)
            logger.info(f"Unregistered calculation {calculation_id} - freed {calculation.cpu_cores} CPU cores, "
                       f"{calculation.estimated_memory_mb} MB memory")
    
    def get_active_calculations(self) -> List[CalculationResourceUsage]:
        """Get list of active calculations and their resource usage."""
        return list(self._active_calculations.values())
    
    def get_total_allocated_resources(self) -> Tuple[int, int]:
        """Get total allocated CPU cores and memory."""
        total_cpu_cores = sum(calc.cpu_cores for calc in self._active_calculations.values())
        total_memory_mb = sum(calc.estimated_memory_mb for calc in self._active_calculations.values())
        return total_cpu_cores, total_memory_mb
    
    def can_allocate_resources(self, cpu_cores: int, memory_mb: int, calculation_method: str = "DFT") -> Tuple[bool, str]:
        """
        Check if resources can be allocated for a new calculation.
        
        Returns:
            Tuple of (can_allocate: bool, reason: str)
        """
        estimated_memory = self._estimate_memory_usage(memory_mb, calculation_method, cpu_cores)
        
        # Get current system state
        system_info = self.get_system_info()
        allocated_cpu, allocated_memory = self.get_total_allocated_resources()
        
        # Calculate resource usage after allocation
        total_cpu_after = allocated_cpu + cpu_cores
        total_memory_after = allocated_memory + estimated_memory
        
        # Check CPU constraints
        max_allowed_cpu = int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0)
        if total_cpu_after > max_allowed_cpu:
            return False, f"CPU cores limit exceeded. Requested: {cpu_cores}, Available: {max_allowed_cpu - allocated_cpu}, Current usage: {allocated_cpu}/{max_allowed_cpu}"
        
        # Check memory constraints
        max_allowed_memory = int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0)
        if total_memory_after > max_allowed_memory:
            return False, f"Memory limit exceeded. Requested: {estimated_memory} MB, Available: {max_allowed_memory - allocated_memory} MB, Current usage: {allocated_memory}/{max_allowed_memory} MB"
        
        # Check current system load (if psutil is available)
        if psutil is not None:
            if system_info.cpu_usage_percent > self._resource_constraints.max_cpu_utilization_percent:
                return False, f"System CPU usage too high: {system_info.cpu_usage_percent:.1f}% > {self._resource_constraints.max_cpu_utilization_percent}%"
            
            if system_info.memory_usage_percent > self._resource_constraints.max_memory_utilization_percent:
                return False, f"System memory usage too high: {system_info.memory_usage_percent:.1f}% > {self._resource_constraints.max_memory_utilization_percent}%"
        
        return True, "Resources available"
    
    def get_resource_summary(self) -> Dict:
        """Get a summary of resource usage and constraints."""
        system_info = self.get_system_info()
        allocated_cpu, allocated_memory = self.get_total_allocated_resources()
        
        max_allowed_cpu = int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0)
        max_allowed_memory = int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0)
        
        return {
            "system_info": {
                "total_cpu_cores": system_info.total_cpu_cores,
                "total_memory_mb": system_info.total_memory_mb,
                "available_memory_mb": system_info.available_memory_mb,
                "cpu_usage_percent": system_info.cpu_usage_percent,
                "memory_usage_percent": system_info.memory_usage_percent,
                "timestamp": system_info.timestamp.isoformat()
            },
            "resource_constraints": {
                "max_cpu_utilization_percent": self._resource_constraints.max_cpu_utilization_percent,
                "max_memory_utilization_percent": self._resource_constraints.max_memory_utilization_percent,
                "max_allowed_cpu_cores": max_allowed_cpu,
                "max_allowed_memory_mb": max_allowed_memory
            },
            "allocated_resources": {
                "total_allocated_cpu_cores": allocated_cpu,
                "total_allocated_memory_mb": allocated_memory,
                "available_cpu_cores": max_allowed_cpu - allocated_cpu,
                "available_memory_mb": max_allowed_memory - allocated_memory,
                "active_calculations_count": len(self._active_calculations)
            }
        }
    
    def _estimate_memory_usage(self, requested_memory_mb: int, calculation_method: str, cpu_cores: int) -> int:
        """
        Estimate actual memory usage based on calculation method and parameters.
        
        This accounts for overhead, parallel processing, and method-specific requirements.
        """
        base_memory = max(requested_memory_mb, 512)  # Minimum 512 MB
        
        # Method-specific multipliers
        method_multipliers = {
            "DFT": 1.2,      # 20% overhead for DFT calculations
            "HF": 1.1,       # 10% overhead for HF calculations
            "MP2": 1.5,      # 50% overhead for MP2 (higher memory usage)
            "CCSD": 2.0,     # 100% overhead for CCSD (very memory intensive)
            "CCSD_T": 2.5,   # 150% overhead for CCSD(T)
            "TDDFT": 1.4     # 40% overhead for TDDFT
        }
        
        multiplier = method_multipliers.get(calculation_method, 1.2)
        
        # Additional overhead for parallel processing
        if cpu_cores > 1:
            # Add 10% per additional core (communication overhead)
            parallel_overhead = 1.0 + 0.1 * (cpu_cores - 1)
            multiplier *= parallel_overhead
        
        # System overhead (temporary files, OS buffers, etc.)
        system_overhead = 1.1
        
        estimated_memory = int(base_memory * multiplier * system_overhead)
        
        # Cap estimate at reasonable maximum (don't go over 4x requested)
        max_estimate = base_memory * 4
        estimated_memory = min(estimated_memory, max_estimate)
        
        return estimated_memory


# Global resource manager instance
_resource_manager: Optional[SystemResourceManager] = None


def get_resource_manager() -> SystemResourceManager:
    """Get the global resource manager instance."""
    global _resource_manager
    if _resource_manager is None:
        _resource_manager = SystemResourceManager()
    return _resource_manager


def shutdown_resource_manager():
    """Shutdown the global resource manager."""
    global _resource_manager
    if _resource_manager is not None:
        logger.info("Shutting down resource manager")
        _resource_manager = None