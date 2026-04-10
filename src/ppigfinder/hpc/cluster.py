"""
HPC Cluster Management Module.
Handles SSH connections, SFTP file transfers, and Slurm job submissions.
"""
import paramiko
import logging
from typing import Optional

logger = logging.getLogger(__name__)

class ClusterManager:
    def __init__(self, host: str, username: str, key_filepath: Optional[str] = None):
        """Initialises the SSH client configuration."""
        self.host = host
        self.username = username
        self.key_filepath = key_filepath
        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
    def connect(self) -> None:
        """Establishes an SSH connection to the cluster."""
        try:
            if self.key_filepath:
                self.client.connect(self.host, username=self.username, key_filename=self.key_filepath)
            else:
                # Relies on ~/.ssh/config or standard key locations
                self.client.connect(self.host, username=self.username)
            logger.info(f"Successfully connected to {self.host}")
        except Exception as e:
            logger.error(f"Failed to connect to {self.host}: {e}")
            raise

    def upload_file(self, local_path: str, remote_path: str) -> None:
        """Uploads a file to the cluster via SFTP."""
        sftp = self.client.open_sftp()
        try:
            sftp.put(local_path, remote_path)
            logger.info(f"Uploaded {local_path} to {remote_path}")
        finally:
            sftp.close()

    def submit_slurm_job(self, sbatch_script_path: str) -> str:
        """Submits a job to Slurm and returns the Job ID."""
        command = f"sbatch {sbatch_script_path}"
        stdin, stdout, stderr = self.client.exec_command(command)
        
        output = stdout.read().decode().strip()
        error = stderr.read().decode().strip()
        
        if error:
            logger.error(f"Slurm submission error: {error}")
            raise RuntimeError(error)
            
        return output

    def close(self) -> None:
        """Closes the SSH connection."""
        self.client.close()
