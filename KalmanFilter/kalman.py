from env import *
class KalmanFilter:
    def __init__(self, noise_velocity, noise_position) -> None:
        # Complete this function to construct

        # Assume that nothing is known 
        # about the state of the target at this instance
        
        pass

    def input(self, observed_state:State, accel:numpy.ndarray, justUpdated:bool):

        # This function is executed multiple times during the reading.
        # When an observation is read, the `justUpdated` is true, otherwise it is false
        
        # accel is the acceleration(control) vector. 
        # It is dynamically obtained regardless of the state of the RADAR 
        # (i.e regardless of `justUpdated`) 

        # When `justUpdated` is false, the state is the same as the previously provided state
        # (i.e, given `observed_state` is not updated, it's the same as the previous outdated one)


        # Complete this function where current estimate of target is updated
        pass

    def get_current_estimate(self)->State:
        
        # Complete this function where the current state of the target is returned
        
        pass