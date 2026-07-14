import os
import numpy as np
import pandas as pd

class VelocityDataLoader:
    """
    Responsible for loading and cleaning the input velocity field data
    from the given file.
    """

    def __init__(self, filepath: str):
        self.filepath = filepath

    def load_data(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Reads the data and extracts the necessary arrays.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: 
                x, y, u, v arrays extracted from the dataset.
                
        Raises:
            FileNotFoundError: If the data file does not exist.
        """
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f"Velocity data file not found: {self.filepath}")

        df = pd.read_csv(self.filepath)
        df.columns = df.columns.str.strip()

        x = df['x-coordinate'].to_numpy()
        y = df['y-coordinate'].to_numpy()
        u = df['x-velocity'].to_numpy()
        v = df['y-velocity'].to_numpy()

        return x, y, u, v
