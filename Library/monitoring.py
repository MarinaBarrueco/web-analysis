"""
Runtime Monitoring and Safety System for Peptide Analysis Pipeline
"""

import time
import psutil
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional
from pathlib import Path
import warnings
from contextlib import contextmanager

class AnalysisMonitor:
    """Comprehensive monitoring system for the analysis pipeline"""
    
    def __init__(self, log_file: Optional[str] = None):
        self.start_time = time.time()
        self.step_times = {}
        self.memory_usage = []
        self.warnings_count = 0
        self.errors_count = 0
        self.data_quality_metrics = {}
        
        # Setup logging
        if log_file:
            logging.basicConfig(
                filename=log_file,
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s'
            )
        else:
            logging.basicConfig(level=logging.INFO)
        
        self.logger = logging.getLogger(__name__)
        self.logger.info("Analysis monitor initialized")
    
    @contextmanager
    def monitor_step(self, step_name: str):
        """Context manager to monitor individual analysis steps"""
        start_time = time.time()
        start_memory = psutil.virtual_memory().used / 1024 / 1024  # MB
        
        self.logger.info(f"Starting step: {step_name}")
        
        try:
            yield self
        except Exception as e:
            self.errors_count += 1
            self.logger.error(f"Error in step {step_name}: {str(e)}")
            raise
        finally:
            end_time = time.time()
            end_memory = psutil.virtual_memory().used / 1024 / 1024  # MB
            
            duration = end_time - start_time
            memory_delta = end_memory - start_memory
            
            self.step_times[step_name] = duration
            self.memory_usage.append({
                'step': step_name,
                'duration': duration,
                'memory_used_mb': end_memory,
                'memory_delta_mb': memory_delta
            })
            
            self.logger.info(
                f"Completed step: {step_name} "
                f"(Duration: {duration:.2f}s, Memory: +{memory_delta:.1f}MB)"
            )
    
    def check_data_quality(self, df: pd.DataFrame, step_name: str) -> Dict[str, Any]:
        """Perform comprehensive data quality checks"""
        quality_metrics = {
            'step': step_name,
            'timestamp': time.time(),
            'rows': len(df),
            'columns': len(df.columns),
            'memory_usage_mb': df.memory_usage(deep=True).sum() / 1024 / 1024,
            'missing_values': df.isnull().sum().sum(),
            'duplicate_rows': df.duplicated().sum(),
            'warnings': []
        }
        
        # Check for suspiciously small datasets
        if len(df) < 10:
            warning = f"Very small dataset in {step_name}: {len(df)} rows"
            quality_metrics['warnings'].append(warning)
            self.logger.warning(warning)
            self.warnings_count += 1
        
        # Check for high missing value percentage
        missing_pct = (quality_metrics['missing_values'] / (len(df) * len(df.columns))) * 100
        if missing_pct > 50:
            warning = f"High missing values in {step_name}: {missing_pct:.1f}%"
            quality_metrics['warnings'].append(warning)
            self.logger.warning(warning)
            self.warnings_count += 1
        
        # Check for numeric columns with suspicious values
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            if (df[col] < 0).any():
                warning = f"Negative values found in {col} at step {step_name}"
                quality_metrics['warnings'].append(warning)
                self.logger.warning(warning)
                self.warnings_count += 1
            
            if df[col].std() == 0:
                warning = f"Zero variance in column {col} at step {step_name}"
                quality_metrics['warnings'].append(warning)
                self.logger.warning(warning)
                self.warnings_count += 1
        
        self.data_quality_metrics[step_name] = quality_metrics
        return quality_metrics
    
    def check_system_resources(self) -> Dict[str, Any]:
        """Monitor system resource usage"""
        cpu_percent = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')
        
        resources = {
            'cpu_percent': cpu_percent,
            'memory_total_gb': memory.total / 1024 / 1024 / 1024,
            'memory_used_gb': memory.used / 1024 / 1024 / 1024,
            'memory_percent': memory.percent,
            'disk_free_gb': disk.free / 1024 / 1024 / 1024,
            'disk_percent': (disk.used / disk.total) * 100
        }
        
        # Check for resource warnings
        if memory.percent > 90:
            warning = f"High memory usage: {memory.percent:.1f}%"
            self.logger.warning(warning)
            self.warnings_count += 1
        
        if cpu_percent > 95:
            warning = f"High CPU usage: {cpu_percent:.1f}%"
            self.logger.warning(warning)
            self.warnings_count += 1
        
        if resources['disk_percent'] > 95:
            warning = f"Low disk space: {resources['disk_percent']:.1f}% used"
            self.logger.warning(warning)
            self.warnings_count += 1
        
        return resources
    
    def validate_statistical_assumptions(self, data: pd.DataFrame, 
                                       condition_cols: List[str]) -> List[str]:
        """Check statistical assumptions for differential expression analysis"""
        warnings = []
        
        # Check for adequate sample sizes
        for col in condition_cols:
            if col in data.columns:
                non_zero_count = (data[col] > 0).sum()
                if non_zero_count < len(data) * 0.1:  # Less than 10% non-zero
                    warnings.append(f"Low detection rate in {col}: {non_zero_count}/{len(data)}")
        
        # Check for sufficient dynamic range
        for col in condition_cols:
            if col in data.columns:
                col_data = data[col]
                if col_data.max() / (col_data.min() + 1) < 10:  # Less than 10-fold range
                    warnings.append(f"Limited dynamic range in {col}")
        
        # Check for outliers
        for col in condition_cols:
            if col in data.columns:
                col_data = data[col]
                q75, q25 = np.percentile(col_data, [75, 25])
                iqr = q75 - q25
                outliers = ((col_data < (q25 - 1.5 * iqr)) | 
                           (col_data > (q75 + 1.5 * iqr))).sum()
                if outliers > len(data) * 0.05:  # More than 5% outliers
                    warnings.append(f"High outlier rate in {col}: {outliers}/{len(data)}")
        
        if warnings:
            self.warnings_count += len(warnings)
            for warning in warnings:
                self.logger.warning(f"Statistical assumption warning: {warning}")
        
        return warnings
    
    def generate_safety_report(self) -> Dict[str, Any]:
        """Generate comprehensive safety and performance report"""
        total_time = time.time() - self.start_time
        
        report = {
            'analysis_summary': {
                'total_duration_seconds': total_time,
                'total_steps': len(self.step_times),
                'warnings_count': self.warnings_count,
                'errors_count': self.errors_count,
                'status': 'COMPLETED' if self.errors_count == 0 else 'COMPLETED_WITH_ERRORS'
            },
            'performance_metrics': {
                'step_times': self.step_times,
                'memory_usage': self.memory_usage,
                'slowest_step': max(self.step_times.items(), key=lambda x: x[1]) if self.step_times else None,
                'peak_memory_mb': max([m['memory_used_mb'] for m in self.memory_usage]) if self.memory_usage else 0
            },
            'data_quality': self.data_quality_metrics,
            'system_resources': self.check_system_resources(),
            'recommendations': self._generate_recommendations()
        }
        
        self.logger.info(f"Analysis completed in {total_time:.2f} seconds with {self.warnings_count} warnings")
        return report
    
    def _generate_recommendations(self) -> List[str]:
        """Generate recommendations based on monitoring data"""
        recommendations = []
        
        # Performance recommendations
        if self.step_times:
            longest_step = max(self.step_times.items(), key=lambda x: x[1])
            if longest_step[1] > 300:  # More than 5 minutes
                recommendations.append(
                    f"Consider optimizing '{longest_step[0]}' step which took {longest_step[1]:.1f} seconds"
                )
        
        # Memory recommendations
        if self.memory_usage:
            peak_memory = max([m['memory_used_mb'] for m in self.memory_usage])
            if peak_memory > 8000:  # More than 8GB
                recommendations.append(
                    f"High memory usage detected ({peak_memory:.0f}MB). Consider processing smaller batches."
                )
        
        # Data quality recommendations  
        if self.warnings_count > 10:
            recommendations.append(
                f"High number of warnings ({self.warnings_count}). Review data quality and preprocessing steps."
            )
        
        # General recommendations
        if not recommendations:
            recommendations.append("Analysis completed successfully with no major issues detected.")
        
        return recommendations
    
    def save_report(self, output_path: str):
        """Save monitoring report to file"""
        report = self.generate_safety_report()
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        import json
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        self.logger.info(f"Safety report saved to {output_path}")

# Utility functions for streamlit integration
def create_streamlit_monitor():
    """Create a monitor instance for Streamlit applications"""
    return AnalysisMonitor()

def display_monitoring_info(monitor: AnalysisMonitor, streamlit_obj):
    """Display monitoring information in Streamlit"""
    if monitor.warnings_count > 0:
        streamlit_obj.warning(f"⚠️ {monitor.warnings_count} warnings detected during analysis")
    
    if monitor.errors_count > 0:
        streamlit_obj.error(f"❌ {monitor.errors_count} errors encountered")
    
    # Show performance metrics
    if monitor.step_times:
        total_time = sum(monitor.step_times.values())
        streamlit_obj.info(f"⏱️ Analysis completed in {total_time:.2f} seconds")
        
        # Show step breakdown
        if len(monitor.step_times) > 1:
            step_df = pd.DataFrame([
                {'Step': step, 'Duration (s)': duration}
                for step, duration in monitor.step_times.items()
            ])
            streamlit_obj.dataframe(step_df)

# Context manager for safe analysis execution
@contextmanager
def safe_analysis_execution(streamlit_obj=None):
    """Context manager for safe analysis execution with monitoring"""
    monitor = create_streamlit_monitor()
    
    try:
        yield monitor
    except Exception as e:
        monitor.errors_count += 1
        monitor.logger.error(f"Analysis failed: {str(e)}")
        if streamlit_obj:
            streamlit_obj.error(f"Analysis failed: {str(e)}")
        raise
    finally:
        # Display final monitoring information
        if streamlit_obj:
            display_monitoring_info(monitor, streamlit_obj)
        
        # Generate and log final report
        report = monitor.generate_safety_report()
        monitor.logger.info("Analysis monitoring completed")