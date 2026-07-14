from dataclasses import dataclass

@dataclass
class BoundaryConditions:
    """
    Encapsulates the Dirichlet boundary conditions for a scalar.
    """
    west: float
    east: float
    south: float
    north: float

    @classmethod
    def from_dict(cls, bc_dict: dict) -> 'BoundaryConditions':
        """
        Creates an instance from a dictionary.
        """
        return cls(
            west=bc_dict.get('west', 0.0),
            east=bc_dict.get('east', 0.0),
            south=bc_dict.get('south', 0.0),
            north=bc_dict.get('north', 0.0)
        )
