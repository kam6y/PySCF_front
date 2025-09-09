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
        self._last_system_info: Optional[SystemResourceInfo] = None
        self._psutil_available = psutil is not None
        self._fallback_mode = False
        
        # Initialize system information with fallback handling
        try:
            if psutil is not None:
                # Get system information using psutil
                self._resource_constraints.system_total_cores = psutil.cpu_count(logical=True)
                self._resource_constraints.system_total_memory_mb = int(psutil.virtual_memory().total / (1024 * 1024))
                logger.info("SystemResourceManager initialized with psutil")
            else:
                raise ImportError("psutil not available")
        except (ImportError, Exception) as e:
            logger.warning(f"Failed to initialize with psutil: {e}. Using fallback mode.")
            self._fallback_mode = True
            self._psutil_available = False
            
            # Fallback to basic multiprocessing info
            try:
                self._resource_constraints.system_total_cores = multiprocessing.cpu_count()
            except Exception:
                # Ultimate fallback
                self._resource_constraints.system_total_cores = 4
                logger.warning("Failed to detect CPU count, using fallback value of 4 cores")
            
            # Conservative memory estimate
            self._resource_constraints.system_total_memory_mb = 4096
            logger.warning("Using conservative memory estimate of 4GB")
            
        logger.info(f"Initialized SystemResourceManager - Mode: {'psutil' if self._psutil_available else 'fallback'}, "
                   f"CPU cores: {self._resource_constraints.system_total_cores}, "
                   f"Total memory: {self._resource_constraints.system_total_memory_mb} MB")
    
    def get_system_info(self) -> SystemResourceInfo:
        """Get current system resource information."""
        current_time = datetime.now(timezone.utc)
        
        if self._fallback_mode or not self._psutil_available:
            # Return conservative estimates when psutil is not available
            system_info = SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=self._resource_constraints.system_total_memory_mb,
                available_memory_mb=int(self._resource_constraints.system_total_memory_mb * 0.6),  # Assume 60% available
                cpu_usage_percent=30.0,  # Conservative estimate
                memory_usage_percent=40.0,  # Conservative estimate
                timestamp=current_time
            )
            logger.debug("Using fallback system info")
            return system_info
        
        # Try to get system info using psutil with fallback handling
        try:
            memory = psutil.virtual_memory()
            cpu_percent = psutil.cpu_percent(interval=0.1)
            
            system_info = SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=int(memory.total / (1024 * 1024)),
                available_memory_mb=int(memory.available / (1024 * 1024)),
                cpu_usage_percent=cpu_percent,
                memory_usage_percent=memory.percent,
                timestamp=current_time
            )
            logger.debug(f"Got system info via psutil: CPU {cpu_percent:.1f}%, Memory {memory.percent:.1f}%")
            return system_info
            
        except Exception as e:
            logger.warning(f"Error getting system information via psutil: {e}. Switching to fallback mode.")
            self._fallback_mode = True
            self._psutil_available = False
            
            # Return fallback values
            system_info = SystemResourceInfo(
                total_cpu_cores=self._resource_constraints.system_total_cores,
                total_memory_mb=self._resource_constraints.system_total_memory_mb,
                available_memory_mb=int(self._resource_constraints.system_total_memory_mb * 0.6),
                cpu_usage_percent=30.0,  # Conservative estimate
                memory_usage_percent=40.0,  # Conservative estimate  
                timestamp=current_time
            )
            logger.debug("Using fallback system info due to psutil error")
            return system_info
    
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
        
        # Get current resource state for logging
        allocated_cpu, allocated_memory = self.get_total_allocated_resources()
        
        self._active_calculations[calculation_id] = calculation_usage
        
        # Log detailed resource allocation
        new_allocated_cpu = allocated_cpu + cpu_cores
        new_allocated_memory = allocated_memory + estimated_memory
        logger.info(f"Registered calculation {calculation_id} - CPU cores: {cpu_cores}, "
                   f"Memory: {memory_mb} MB (estimated: {estimated_memory} MB)")
        logger.info(f"Total allocated resources after registration: CPU {new_allocated_cpu}/{self._resource_constraints.system_total_cores}, "
                   f"Memory {new_allocated_memory}/{self._resource_constraints.system_total_memory_mb} MB")
    
    def unregister_calculation(self, calculation_id: str) -> None:
        """Unregister a calculation when it completes."""
        if calculation_id in self._active_calculations:
            calculation = self._active_calculations.pop(calculation_id)
            
            # Get remaining resource state for logging
            remaining_cpu, remaining_memory = self.get_total_allocated_resources()
            
            logger.info(f"Unregistered calculation {calculation_id} - freed {calculation.cpu_cores} CPU cores, "
                       f"{calculation.estimated_memory_mb} MB memory")
            logger.info(f"Remaining allocated resources: CPU {remaining_cpu}/{self._resource_constraints.system_total_cores}, "
                       f"Memory {remaining_memory}/{self._resource_constraints.system_total_memory_mb} MB")
        else:
            logger.warning(f"Attempted to unregister calculation {calculation_id} that was not registered")
    
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
        
        # Log detailed resource check
        logger.debug(f"Resource allocation check for {calculation_method}: "
                    f"Requested CPU: {cpu_cores}, Memory: {memory_mb} MB (estimated: {estimated_memory} MB)")
        logger.debug(f"Current allocation: CPU {allocated_cpu}, Memory {allocated_memory} MB, Active calculations: {len(self._active_calculations)}")
        
        # Check CPU constraints
        max_allowed_cpu = max(1, int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0))
        if total_cpu_after > max_allowed_cpu:
            reason = f"CPU cores limit exceeded. Requested: {cpu_cores}, Available: {max_allowed_cpu - allocated_cpu}, Current usage: {allocated_cpu}/{max_allowed_cpu}"
            logger.debug(f"Resource allocation failed: {reason}")
            return False, reason
        
        # Check memory constraints
        max_allowed_memory = max(256, int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0))
        if total_memory_after > max_allowed_memory:
            reason = f"Memory limit exceeded. Requested: {estimated_memory} MB, Available: {max_allowed_memory - allocated_memory} MB, Current usage: {allocated_memory}/{max_allowed_memory} MB"
            logger.debug(f"Resource allocation failed: {reason}")
            return False, reason
        
        # Check current system load (if psutil is available and not in fallback mode)
        if not self._fallback_mode and self._psutil_available:
            try:
                if system_info.cpu_usage_percent > self._resource_constraints.max_cpu_utilization_percent:
                    reason = f"System CPU usage too high: {system_info.cpu_usage_percent:.1f}% > {self._resource_constraints.max_cpu_utilization_percent}%"
                    logger.debug(f"Resource allocation failed: {reason}")
                    return False, reason
                
                if system_info.memory_usage_percent > self._resource_constraints.max_memory_utilization_percent:
                    reason = f"System memory usage too high: {system_info.memory_usage_percent:.1f}% > {self._resource_constraints.max_memory_utilization_percent}%"
                    logger.debug(f"Resource allocation failed: {reason}")
                    return False, reason
            except Exception as e:
                logger.warning(f"Failed to check system load: {e}. Proceeding with resource allocation.")
        else:
            logger.debug("Skipping system load check - fallback mode active")
        
        success_reason = f"Resources available - CPU: {total_cpu_after}/{max_allowed_cpu}, Memory: {total_memory_after}/{max_allowed_memory} MB"
        logger.debug(f"Resource allocation succeeded: {success_reason}")
        return True, success_reason
    
    def get_resource_summary(self) -> Dict:
        """Get a summary of resource usage and constraints."""
        system_info = self.get_system_info()
        allocated_cpu, allocated_memory = self.get_total_allocated_resources()
        
        max_allowed_cpu = max(1, int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0))
        max_allowed_memory = max(256, int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0))
        
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
    
    def has_resources_improved(self, improvement_threshold_percent: float = 10.0) -> Tuple[bool, str]:
        """
        Check if system resources have improved significantly since the last check.
        
        Args:
            improvement_threshold_percent: Minimum improvement threshold (default: 10%)
            
        Returns:
            Tuple of (has_improved: bool, reason: str)
        """
        if self._fallback_mode or not self._psutil_available:
            # In fallback mode, periodically return True to allow queue processing
            # This prevents calculations from being stuck indefinitely
            logger.debug("Resource improvement check skipped - fallback mode active")
            return True, "Fallback mode - allowing queue processing"
        
        try:
            if self._last_system_info is None:
                # No previous data to compare - consider improved
                return True, "No previous resource data to compare"
            
            current_info = self.get_system_info()
            previous_info = self._last_system_info
            
            # Check CPU usage improvement
            cpu_improvement = previous_info.cpu_usage_percent - current_info.cpu_usage_percent
            
            # Check memory usage improvement  
            memory_improvement = previous_info.memory_usage_percent - current_info.memory_usage_percent
            
            # Consider improved if either CPU or memory usage dropped by threshold amount
            cpu_improved = cpu_improvement >= improvement_threshold_percent
            memory_improved = memory_improvement >= improvement_threshold_percent
            
            if cpu_improved or memory_improved:
                improvements = []
                if cpu_improved:
                    improvements.append(f"CPU usage decreased by {cpu_improvement:.1f}%")
                if memory_improved:
                    improvements.append(f"memory usage decreased by {memory_improvement:.1f}%")
                
                reason = f"Resources improved: {', '.join(improvements)}"
                logger.info(reason)
                return True, reason
            else:
                reason = f"No significant resource improvement (CPU: {cpu_improvement:+.1f}%, Memory: {memory_improvement:+.1f}%)"
                return False, reason
                
        except Exception as e:
            logger.warning(f"Failed to check resource improvements: {e}")
            # In case of error, allow queue processing to prevent calculations from being stuck
            return True, f"Resource improvement check failed - allowing queue processing"
    
    def check_if_system_resources_insufficient(self) -> Tuple[bool, str]:
        """
        Check if system resources are fundamentally insufficient even with no running calculations.
        
        Returns:
            Tuple of (insufficient: bool, reason: str)
        """
        if self._fallback_mode or not self._psutil_available:
            # In fallback mode, assume resources are sufficient unless overridden
            # This prevents false errors when system monitoring is not available
            logger.debug("Resource insufficiency check skipped - fallback mode active")
            return False, "Resource monitoring unavailable - assuming sufficient resources"
            
        try:
            system_info = self.get_system_info()
            
            # Check if current system load exceeds limits even with no active calculations
            if len(self._active_calculations) == 0:
                if system_info.cpu_usage_percent > self._resource_constraints.max_cpu_utilization_percent:
                    return True, f"System CPU usage ({system_info.cpu_usage_percent:.1f}%) exceeds limit ({self._resource_constraints.max_cpu_utilization_percent}%) with no active calculations"
                    
                if system_info.memory_usage_percent > self._resource_constraints.max_memory_utilization_percent:
                    return True, f"System memory usage ({system_info.memory_usage_percent:.1f}%) exceeds limit ({self._resource_constraints.max_memory_utilization_percent}%) with no active calculations"
            
            return False, "System resources are sufficient"
            
        except Exception as e:
            logger.warning(f"Failed to check system resource insufficiency: {e}")
            # In case of error, assume resources are sufficient to avoid blocking calculations
            return False, f"Resource check failed ({str(e)}) - assuming sufficient resources"
    
    def get_diagnostics(self) -> Dict:
        """Get diagnostic information about the resource manager state."""
        try:
            system_info = self.get_system_info()
            allocated_cpu, allocated_memory = self.get_total_allocated_resources()
            
            return {
                "mode": "psutil" if (self._psutil_available and not self._fallback_mode) else "fallback",
                "psutil_available": self._psutil_available,
                "fallback_mode": self._fallback_mode,
                "system_info": {
                    "total_cpu_cores": system_info.total_cpu_cores,
                    "total_memory_mb": system_info.total_memory_mb,
                    "current_cpu_usage_percent": system_info.cpu_usage_percent,
                    "current_memory_usage_percent": system_info.memory_usage_percent,
                    "available_memory_mb": system_info.available_memory_mb,
                },
                "resource_constraints": {
                    "max_cpu_utilization_percent": self._resource_constraints.max_cpu_utilization_percent,
                    "max_memory_utilization_percent": self._resource_constraints.max_memory_utilization_percent,
                },
                "active_calculations": {
                    "count": len(self._active_calculations),
                    "allocated_cpu_cores": allocated_cpu,
                    "allocated_memory_mb": allocated_memory,
                    "calculation_details": [
                        {
                            "id": calc.calculation_id,
                            "cpu_cores": calc.cpu_cores,
                            "memory_mb": calc.memory_mb,
                            "estimated_memory_mb": calc.estimated_memory_mb
                        }
                        for calc in self._active_calculations.values()
                    ]
                },
                "capacity": {
                    "max_allowed_cpu": max(1, int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0)),
                    "max_allowed_memory_mb": max(256, int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0)),
                    "available_cpu_cores": max(0, max(1, int(system_info.total_cpu_cores * self._resource_constraints.max_cpu_utilization_percent / 100.0)) - allocated_cpu),
                    "available_memory_mb": max(0, max(256, int(system_info.total_memory_mb * self._resource_constraints.max_memory_utilization_percent / 100.0)) - allocated_memory)
                }
            }
        except Exception as e:
            logger.error(f"Failed to get resource diagnostics: {e}")
            return {
                "error": str(e),
                "mode": "error",
                "psutil_available": self._psutil_available,
                "fallback_mode": self._fallback_mode
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