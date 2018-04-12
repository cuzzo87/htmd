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

        if initReps and isinstance(viewer, VMD):
            viewer.send('color Display Background {}'.format(self.background))
            viewer.loadMol(mol)
            mol.reps._repsVMD(viewer)

    def saveScene(self, repr='current'):
        # repr current, previous
        viewer = self.viewer
        mol = self.mol

        if isinstance(viewer, VMD):
            rotate_matrix = self._getMatrixVMD('rotate')
            center_matrix = self._getMatrixVMD('center')
            scale_matrix = self._getMatrixVMD('scale')
            global_matrix = self._getMatrixVMD('global')

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

        viewer = self.viewer

        rot_mat = str(scene.rotate_matrix.tolist())
        rot_mat = rot_mat.replace('[', '{').replace(']', '}').replace(',', '')

        viewer.send('molinfo top set rotate_matrix [list {}]'.format(rot_mat))


    def createAnimation(self, startSceneId, endSceneId):

        start_scene = self.scenes[startSceneId]
        start_rotMat = start_scene.rotate_matrix




        end_scene = self.scenes[endSceneId]
        end_rotMat = end_scene.rotate_matrix


    def _getMatrixVMD(self, matrixtype):
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

def _bestEuler(euler):
    from math import pi

    new_euler = []

    for e in euler:
        if e > pi:
            new_euler.append(e-2*pi)
        else:
            new_euler.append(e + 2 * pi)
    return np.array(new_euler)
