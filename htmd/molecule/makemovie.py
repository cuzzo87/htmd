from htmd.vmdviewer import VMD, getCurrentViewer
from nglview import NGLWidget
from tempfile import NamedTemporaryFile
import numpy as np

class Scene:

    def __init__(self, reprs=None, rot_mat=None):

        self.representations = reprs
        self.rotate_matrix = rot_mat

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
            rotate_matrix = self._getRotateMatrixVMD()

        #TODO ngl

        scene = Scene(reprs=mol.reps, rot_mat=rotate_matrix)

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

        s





    def _getRotateMatrixVMD(self):
        viewer = self.viewer
        outputFile = NamedTemporaryFile(delete=False).name

        print('save scene for vmd ')
        viewer.send('set R [molinfo top get rotate_matrix]')


        self._writeTclOutput('R', outputFile)

        f = open(outputFile, 'r')
        txt_rotateMatrix = f.read().replace('{', '').replace('}', '')
        rotateMatrix = np.array(txt_rotateMatrix.split(), dtype=float)
        rotateMatrix = rotateMatrix.reshape(4,4)

        return rotateMatrix

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

