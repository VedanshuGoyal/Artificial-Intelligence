from env import *
import numpy as np

class KalmanFilter:

    def __init__(self, noise_velocity, noise_position) -> None:
        # Complete this function to construct

        # Assume that nothing is known 
        # about the state of the target at this instance

        # making A and B arrays
        
        _D = self._D = 3 # Dimensions
        _dt = self._dt = 1 # processor clock speed

        self._A = np.eye(2*_D)
        for i in range(0, _D):
            self._A[i][_D + i] = _dt

        self._B = np.zeros((2*_D, _D))
        for i in range(0, _D):
            self._B[i][i] = .5*_dt*_dt
            self._B[_D + i][i] = _dt

        # noise of measurements..
        self._R = np.eye(2*_D)
        for i in range(0, _D):
            self._R[i][i] = noise_position
            self._R[_D + i][_D + i] = noise_velocity

        self._H = np.eye(2*_D)

        self._isinit = 0;

    def __StateToMatrix(self, state : State):
        x = np.zeros((2*self._D, 1))
        for i in range(self._D):
            x[i] = state.position[i]
            x[self._D + i] = state.velocity[i]
        return x

    def __MatrixToState(self, x):
        return State(x[:self._D, :], x[self._D:, :])

    def __firstTimeinput(self, observed_state:State):
        self._x = self.__StateToMatrix(observed_state)
        self._P = np.copy(self._R)
        _isinit = 1;
        return


    def __predict(self, accel:numpy.ndarray):
        # x = A.x + B.a + W
        # P = A.P.At + Q 
        # Assuming 0 variance in Accel...
        # W & Q = 0

        accel = accel.reshape((3, 1))
        
        new_x = self._A.dot(self._x) + self._B.dot(accel)
        new_P = self._A.dot(self._P).dot(self._A.T)

        self._x = new_x
        self._P = new_P

    def __update(self, z : numpy.ndarray, R : numpy.ndarray):
        # z -> measurement ; R -> Noise
        # y = z - Hx
        # S = H.P.Ht + R
        # Kalman Gain -> K = P. Ht . S^-1
        # x = x + K.y
        # P = (I - K.H) * P

        y = z - self._H.dot(self._x)
        S = self._H.dot(self._P).dot(self._H.T) + R
        K = self._P.dot(self._H.T).dot(np.linalg.inv(S))

        new_x = self._x + K.dot(y)
        new_P = (np.eye(2*self._D) - K.dot(self._H)).dot(self._P)

        self._x = new_x
        self._P = new_P

    def input(self, observed_state:State, accel:numpy.ndarray, justUpdated:bool):

        # This function is executed multiple times during the reading.
        # When an observation is read, the `justUpdated` is true, otherwise it is false

        # accel is the acceleration(control) vector. 
        # It is dynamically obtained regardless of the state of the RADAR 
        # (i.e regardless of `justUpdated`) 

        # When `justUpdated` is false, the state is the same as the previously provided state
        # (i.e, given `observed_state` is not updated, it's the same as the previous outdated one)

        # Complete this function where current estimate of target is updated
        if(self._isinit == 0):
            # first time input
            self.__firstTimeinput(observed_state)
            self._isinit = 1
            return

        self.__predict(accel)

        if(justUpdated == False):
            return

        self.__update(self.__StateToMatrix(observed_state), self._R)


    def get_current_estimate(self)->State:
        
        # Complete this function where the current state of the target is returned
        return self.__MatrixToState(self._x);