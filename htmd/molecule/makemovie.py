from htmd.vmdviewer import VMD, getCurrentViewer
from nglview import NGLWidget
from tempfile import NamedTemporaryFile
import numpy as np

class Scene:

    def __init__(self, reprs=None):

        self.representations = reprs
        self.rotate_matrix = None
        self.center_matrix = None
        self.scale_matrix = None
        self.global_matrix = None

class MovieMaker:

    _stage_parameters = {'background':'white'}

    def __init__(self, mol, viewer=None):

        self.__dict__ = self._stage_parameters
        self.mol = mol
        self.viewer =  self._getViewer(viewer)
        self.scenes = []
        self.animations = []


        self._initStage()

    def _initStage(self):
        mol = self.mol
        viewer = self.viewer

        initReps = False
        if len(mol.reps.replist) != 0:
            initReps = True
        #Todo else

        if initReps and isinstance(viewer, VMD):
            viewer.send('color Display Background {}'.format(self.background))
            viewer.loadMol(mol)
            mol.reps._repsVMD(viewer)

    def saveScene(self, repr='current'):
        # repr current, previous
        viewer = self.viewer
        mol = self.mol

        if isinstance(viewer, VMD):
            rotate_matrix = self._matrixFromVMD('rotate')
            center_matrix = self._matrixFromVMD('center')
            scale_matrix = self._matrixFromVMD('scale')
            global_matrix = self._matrixFromVMD('global')

        #TODO ngl

        scene = Scene(reprs=mol.reps)
        scene.rotate_matrix = rotate_matrix
        scene.center_matrix = center_matrix
        scene.scale_matrix = scale_matrix
        scene.global_matrix = global_matrix

        self.scenes.append(scene)

    def retrieveScene(self, sceneid):
        try:
            scene = self.scenes[sceneid]
        except:
            raise IndexError('The scene with that sceneId does not exists')


        rot_mat = self._matrixToVMD(scene, 'rotate')
        cent_mat = self._matrixToVMD(scene, 'center')
        scale_mat = self._matrixToVMD(scene, 'scale')
        glob_mat = self._matrixToVMD(scene, 'global')


        print(rot_mat)
        viewer = self.viewer

        viewer.send('molinfo top set {{rotate_matrix center_matrix scale_matrix global_matrix}} '
                    '{{{} {} {} {}}}'.format(rot_mat, cent_mat, scale_mat, glob_mat))
        #
        # rot_mat = str(scene.rotate_matrix.tolist())
        # rot_mat = rot_mat.replace('[', '{').replace(']', '}').replace(',', '')
        #
        # viewer.send('molinfo top set rotate_matrix [list {}]'.format(rot_mat))
        # viewer.send('molinfo top set center_matrix [list {}]'.format(rot_mat))
        # viewer.send('molinfo top set scale_matrix [list {}]'.format(rot_mat))
        # viewer.send('molinfo top set global_matrix [list {}]'.format(rot_mat))


    def createAnimation(self, startSceneId, endSceneId):

        start_scene = self.scenes[startSceneId]
        start_rotMat = start_scene.rotate_matrix

        end_scene = self.scenes[endSceneId]
        end_rotMat = end_scene.rotate_matrix

    def _matrixToVMD(self, scene, matrixtype):

        viewer = self.viewer

        mat = getattr(scene, "{}_matrix".format(matrixtype))

        mat = str(mat.tolist())
        mat = mat.replace('[', '{').replace(']', '}').replace(',', '')

        return mat


    def _matrixFromVMD(self, matrixtype):
        viewer = self.viewer
        outputFile = NamedTemporaryFile(delete=False).name

        viewer.send('set R [molinfo top get {}_matrix]'.format(matrixtype))


        self._writeTclOutput('R', outputFile)

        f = open(outputFile, 'r')
        txt_Matrix = f.read().replace('{', '').replace('}', '')
        matrix = np.array(txt_Matrix.split(), dtype=float)
        matrix = matrix.reshape(4,4)

        return matrix


    def _writeTclOutput(self, outputVariable, outputFile):

        viewer = self.viewer

        viewer.send('set Fname [open {} "w"]'.format(outputFile))
        viewer.send('puts -nonewline $Fname ${}'.format(outputVariable) )
        viewer.send('close $Fname')

    def _getViewer(self, viewer):

        if isinstance(viewer, VMD):
            if  viewer.completed():
                return getCurrentViewer()
            return viewer

        elif isinstance(viewer, NGLWidget):
            print('found NGL')
            return viewer

        elif viewer == None:
            return getCurrentViewer()

        else:
            raise ValueError('Not a valid viewer.')




def matrixToEuler(matrix):
    from math import atan2, asin, cos, pi

    m31 = matrix[2][0]
    m12 = matrix[0][1]
    m13 = matrix[0][2]
    m32 = matrix[2][1]
    m33 = matrix[2][2]
    m21 = matrix[1][0]
    m11 = matrix[0][0]

    if m31 == 1:
        phi = 0
        psi = atan2(m12, m13)
        theta = -pi/2
    elif m31 == -1:
        phi = 0
        psi = atan2(m12, m13)
        theta = pi /2

    else:
        theta = -asin(m31)
        cosT = cos(theta)
        psi = atan2(m32/cosT,  m33/cosT)
        phi = atan2(m21/cosT, m11/cosT)

    return np.array([theta, phi, psi])

def matrixToQuaternion(matrix):
    from math import sqrt

    m11 = matrix[0][0]
    m12 = matrix[0][1]
    m13 = matrix[0][2]
    m22 = matrix[1][1]
    m21 = matrix[1][0]
    m23 = matrix[1][2]
    m33 = matrix[2][2]
    m31 = matrix[2][0]
    m32 = matrix[2][1]

    t44as33 = m11 + m22 + m33
    r,s,w,x,y,z = 0,0,0,0,0,0

    if t44as33 > 0:
        r = 1.0 + t44as33
        s = 0.5 / sqrt(r)
        w = s * r
        x = (m23 - m32) * s
        y = (m31 - m13) * s
        z = (m12 - m21) * s

    elif m11 > m22 and m11 > m33:
        r = 1.0 - t44as33 + 2 * m11
        s = 0.5 / sqrt(r)
        w = s * r
        x = (m12 + m21) * s
        y = (m13 + m31) * s
        z = (m23 - m32) * s

    elif m22 > m33:
        r = 1.0 - t44as33 + 2 * m22
        s = 0.5 / sqrt(r)
        w = s * r
        x = (m12 + m21) * s
        y = (m23 + m32) * s
        z = (m31 - m13) * s

    else:
        r = 1.0 - t44as33 + 2 * m33
        s = 0.5 / sqrt(r)
        w = s * r
        x = (m13 + m31) * s
        y = (m23 + m32) * s
        z = (m12 - m21) * s

    return np.array([w,x,y,z])

def quatarc(quat1, quat2, step):
    from math import sqrt, acos, sin

    qdot = np.dot(quat1, quat2)

    if qdot > 0.9999:
        return quat2
    elif qdot < 0:
        quat1 = quat1 * -1

    theta = acos(np.dot(quat1, quat2)/sqrt(np.dot(quat1, quat1) * np.dot(quat2, quat2)))

    quat = np.add( sin(theta * 1-step)/sin(theta) * quat1, sin(theta * step)/sin(theta) * quat2 )

    return quat

def quaternionToMatrix(quat):
    w = quat[0]
    x = quat[1]
    y = quat[2]
    z = quat[3]

    f_row = [1-2*y**2-2*z**2, 2*x*y+2*w*z, 2*x*z-2*w*y, 0]
    s_row = [2*x*y-2*w*z, 1-2*x**2-2*z**2, 2*y*z+2*w*x, 0]
    t_row = [2*x*z+2*w*y, 2*y*z-2*w*x, 1-2*x**2-2*y**2, 0]
    q_row = [0,0,0,1]

    return np.array([f_row, s_row, t_row, q_row])

def _bestEuler(diff_euler, euler):
    from math import pi

    new_euler = []

    for i in range(diff_euler.shape[0]):
        if diff_euler[i] > pi:
            new_euler.append(euler[i]-2*pi)
        else:
            new_euler.append(euler[i] + 2 * pi)
    return np.array(new_euler)
