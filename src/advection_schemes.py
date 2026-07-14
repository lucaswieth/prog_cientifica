from abc import ABC, abstractmethod

class AdvectionScheme(ABC):
    """
    Abstract base class for advection schemes (Strategy Pattern).
    """

    @abstractmethod
    def calculate_coefficients(self, D: float, F: float) -> tuple[float, float]:
        """
        Calculates the matrix coefficients a_neighbor and a_P contributions based on the scheme.

        Args:
            D (float): Diffusion conductance at the face.
            F (float): Convective mass flux at the face.

        Returns:
            tuple[float, float]: (a_neighbor_contribution, a_P_advective_contribution)
        """
        pass

class UpwindScheme(AdvectionScheme):
    """
    Upwind Differencing Scheme (UDS).
    First-order accurate, unconditionally stable, but highly diffusive.
    """

    def calculate_coefficients(self, D: float, F: float, is_positive_direction: bool = True) -> tuple[float, float]:
        """
        Calculates the coefficients using UDS.
        """
        if is_positive_direction:
            # For East/North faces: a_nb uses max(-F, 0), a_P uses + F
            a_nb = D + max(-F, 0.0)
            a_p_contrib = a_nb + F
        else:
            # For West/South faces: a_nb uses max(F, 0), a_P uses - F
            a_nb = D + max(F, 0.0)
            a_p_contrib = a_nb - F
        
        return a_nb, a_p_contrib

class CentralDifferenceScheme(AdvectionScheme):
    """
    Central Differencing Scheme (CDS).
    Second-order accurate, may produce non-physical oscillations if Peclet > 2.
    """

    def calculate_coefficients(self, D: float, F: float, is_positive_direction: bool = True) -> tuple[float, float]:
        """
        Calculates the coefficients using CDS.
        """
        if is_positive_direction:
            a_nb = D - 0.5 * F
            a_p_contrib = a_nb + F
        else:
            a_nb = D + 0.5 * F
            a_p_contrib = a_nb - F
            
        return a_nb, a_p_contrib

def get_advection_scheme(scheme_name: str) -> AdvectionScheme:
    """
    Factory method to instantiate the correct advection scheme.
    """
    schemes = {
        'UDS': UpwindScheme,
        'CDS': CentralDifferenceScheme
    }
    name = scheme_name.upper()
    if name not in schemes:
        raise ValueError(f"Unknown advection scheme '{scheme_name}'. Available: {list(schemes.keys())}")
    
    return schemes[name]()
