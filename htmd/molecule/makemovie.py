from htmd.vmdviewer import VMD, getCurrentViewer
from nglview import NGLWidget,HTMDTrajectory
from tempfile import NamedTemporaryFile
import numpy as np
import os

class Scene:

    def __init__(self, id=None, reprs=None):

        self.id = id

        # vmd
        self.representations = reprs
        self.rotate_matrix = None
        self.center_matrix = None
        self.scale_matrix = None
        self.global_matrix = None

        # ngl
        self.orientation_matrix = None
        self.orientation_delta = None

class Transaction:

    def __init__(self, begScene, endScene, scenes=None):

        self.begScene = begScene
        self.endScene = endScene
        self.begId = self.begScene.id
        self.endId = self.endScene.id

        # vmd
        self.scenes = [] if scenes == None else scenes

        # ngl -->  nglview does not allow to decompose the orientation matrix. Thus we have to store the delta matrix
        #          to apply.



class MovieMaker:

    _stage_parameters = {'background':'white'}

    def __init__(self, mol, viewer=None):

        self.__dict__ = self._stage_parameters
        self.mol = mol
        self.viewer =  self._getViewer(viewer)
        self.scenes = []
        self.animations = []

        self.timeline = []

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
        elif initReps and isinstance(viewer, NGLWidget):
            traj = HTMDTrajectory(mol)
            viewer.add_trajectory(traj)
            mol.reps._repsNGL(viewer)

    def _getScene(self, id=None, repr='current'):
        # repr current, previous
        viewer = self.viewer
        mol = self.mol

        scene = Scene(id=id, reprs=mol.reps)

        if isinstance(viewer, VMD):
            rotate_matrix = self._matrixFromVMD('rotate')
            center_matrix = self._matrixFromVMD('center')
            scale_matrix = self._matrixFromVMD('scale')
            global_matrix = self._matrixFromVMD('global')

            scene.rotate_matrix = rotate_matrix
            scene.center_matrix = center_matrix
            scene.scale_matrix = scale_matrix
            scene.global_matrix = global_matrix

        if isinstance(viewer, NGLWidget):
            orientation_matrix = self._matrixFromNGL()
            scene.orientation_matrix = orientation_matrix

        return scene

    def saveScene(self):
        _id = len(self.scenes)

        scene = self._getScene(id=_id)

        self.scenes.append(scene)

    def retrieveScene(self, scene):
        if isinstance(scene, Scene):
            scene = scene

        else:
            if scene >= len(self.scenes):
                raise IndexError('The scene with that sceneId does not exists')
            scene = self.scenes[scene]

        viewer = self.viewer

        if isinstance(viewer, VMD):
            rot_mat = self._matrixToVMD(scene, 'rotate')
            cent_mat = self._matrixToVMD(scene, 'center')
            scale_mat = self._matrixToVMD(scene, 'scale')
            glob_mat = self._matrixToVMD(scene, 'global')

            _class = viewer
            method = 'send'
            args = 'molinfo top set {{rotate_matrix center_matrix scale_matrix global_matrix}} ' \
                   '{{{} {} {} {}}}'.format(rot_mat, cent_mat, scale_mat, glob_mat)


            self._updateView(_class, method, args)

            #viewer.send('molinfo top set {{rotate_matrix center_matrix scale_matrix global_matrix}} '
             #           '{{{} {} {} {}}}'.format(rot_mat, cent_mat, scale_mat, glob_mat))

        elif isinstance(viewer, NGLWidget):
            _class = viewer.control
            method = 'orient'
            orientation_matrix = self._matrixToNGL(scene)

            self._updateView(_class, method, orientation_matrix)
            viewer.sync_view()

            #viewer._set_camera_orientation(orientation_matrix)
            #viewer.sync_view()

    def _updateView(self, _class, method, args):

        cmd = getattr(_class, method)
        cmd(args)

    def record(self, moviename='mymovie', outdir="imagesVideo", overwrite=False, fps=None, spf=None):
        viewer = self.viewer
        mol = self.mol

        if os.path.isdir(outdir):
            if not overwrite:
                raise IsADirectoryError('The directory exists, change the name or set as True the overwrite argument')
            else:
                os.rmdir(outdir)
        os.mkdir(outdir)


        if isinstance(viewer, VMD):
            self._recordVMD(viewer, mol, outdir)

        elif isinstance(viewer, NGLWidget):
            pass

    def _recordVMD(self, viewer, mol, outdir):

        preload_tcl = """global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
set molID [molinfo top]
set tach "$vmdEnv/tachyon_$Arch"

"""
        viewer.send(preload_tcl)
        frames = mol.numFrames

        for nf in range(frames):
            nframe = "%06d" % nf
            tcl = """animate goto {}
display update ui
set filename {}/image.{}.tga
render Tachyon $filename "$tach" -aasamples 12 %s -format TARGA -o %s
""".format(nf, outdir, nframe)
            viewer.send(tcl)
            if nf == 3:
                break


    def play(self):
        from time import sleep

        viewer = self.viewer

        if isinstance(viewer, NGLWidget):
            print('Not stable. At the moment it does not work')
            return

        for element in self.timeline:
            if isinstance(element, Scene):
                scene = element
                self.retrieveScene(scene)

            if isinstance(element, Transaction):
                transaction = element
                for scene in transaction.scenes:
                    #TODO  temporary approach until nglview allows to decompose the orientation matrix
                    if isinstance(viewer, NGLWidget):
                        _class = viewer.control
                        method = 'apply_matrix'
                        mat = scene.orientation_delta
                        args = np.concatenate(mat).tolist()
                        self._updateView(_class, method, args)
                        viewer.sync_view()
                    else:
                        self.retrieveScene(scene)

    def autoTimeline(self):
        tmp_timeline = []

        transaction_skipped = []

        for i in range(len(self.scenes)):
            scene = self.scenes[i]
            tmp_timeline.append(scene)

            for t in self.animations:
                if t.begId == i  and t not in tmp_timeline:
                    if t.endId == i + 1:
                        tmp_timeline.append(t)
                    else:
                        transaction_skipped.append(t)

        self.timeline = tmp_timeline
        if len(transaction_skipped) != 0:
            print("Some transactions were not possible to place in timeline ( n transaction skipped: {}) ".format(len(transaction_skipped)))




    def transaction(self, begSceneId, endSceneId, numsteps=50):
        #self.retrieveScene(begSceneId)

        begScene = self.scenes[begSceneId]
        endScene = self.scenes[endSceneId]

        viewer = self.viewer

        if isinstance(viewer, VMD):
            listScenes = self._transactionVMD(begScene, endScene, numsteps)
        elif isinstance(viewer, NGLWidget):
            listScenes = self._transactionNGL(begScene, endScene, numsteps)

        T = Transaction(begScene, endScene, listScenes)
        self.animations.append(T)

    def _transactionNGL(self, begScene, endScene, numsteps):

        current_orientation = begScene.orientation_matrix

        beg_orientationQuat = _normalizeQuat(matrixToQuaternion(begScene.orientation_matrix))
        end_orientationQuat = _normalizeQuat(matrixToQuaternion(endScene.orientation_matrix))

        stepsize = 1 / numsteps
        current_deltaQuat = quatarc(beg_orientationQuat, end_orientationQuat, stepsize)

        current_deltaMatrix = quaternionToMatrix(current_deltaQuat)

        list_scenes = []
        for i in range(numsteps+1):
            sc = Scene()
            sc.orientation_delta = current_deltaMatrix

            list_scenes.append(sc)

        return list_scenes

    def _transactionVMD(self, begScene, endScene, numsteps):

        diff_center = endScene.center_matrix - begScene.center_matrix
        diff_scale = endScene.scale_matrix - begScene.scale_matrix
        diff_global = endScene.global_matrix - begScene.global_matrix

        current_rotateMatrix = begScene.rotate_matrix
        current_centerMatrix = begScene.center_matrix
        current_scaleMatrix = begScene.scale_matrix
        current_globalMatrix = begScene.global_matrix

        beg_rotateQuat = matrixToQuaternion(begScene.rotate_matrix)
        end_rotateQuat = matrixToQuaternion(endScene.rotate_matrix)

        stepsize = 1 / numsteps

        list_scenes = []
        for i in range(1, numsteps+1):
            stepratio = stepsize * i
            qarc = quatarc(beg_rotateQuat, end_rotateQuat, stepratio)
            current_rotateMatrix = quaternionToMatrix(qarc)
            current_centerMatrix = np.add(current_centerMatrix, diff_center * stepsize)
            current_scaleMatrix = np.add(current_scaleMatrix, diff_scale * stepsize)
            current_globalMatrix = np.add(current_globalMatrix, diff_global * stepsize)

            sc = Scene()
            sc.rotate_matrix = current_rotateMatrix
            sc.center_matrix = current_centerMatrix
            sc.scale_matrix = current_scaleMatrix
            sc.global_matrix = current_globalMatrix

            list_scenes.append(sc)
        return list_scenes


    def _matrixToVMD(self, scene, matrixtype):

        mat = getattr(scene, "{}_matrix".format(matrixtype))

        mat = str(mat.tolist())
        mat = mat.replace('[', '{').replace(']', '}').replace(',', '')

        return mat

    def _matrixToNGL(self, scene):

        mat = scene.orientation_matrix
        mat = np.concatenate(mat).tolist()

        return mat


    def _matrixFromNGL(self):
        viewer = self.viewer

        rot_matrix = np.array(viewer._camera_orientation).reshape(4,4)
        return rot_matrix


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

def _matrixToQuaternion(matrix):
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
    t = t44as33
    r, s, w, x, y, z = 0, 0, 0, 0, 0, 0

    if t > 0:
        S = sqrt(t+1) * 2

        w = 0.25 * S
        x = (m32 - m32) / S
        y = (m13 - m31) / S
        z = (m21 - m12) / S

    elif m11 > m22 and m11 > m33:
        S = sqrt(1 + m11 - m22 - m33) * 2

        w = (m32 - m23) /S
        x = 0.25 * S
        y = (m12 + m21) /S
        z = (m13 + m31) /S

    elif m22 > m33:
        S = sqrt(1 + m22 - m11 -m33) * 2

        w = (m13 - m31) / S
        x = (m12 + m21) / S
        y = 0.25 * S
        z = (m23 + m32) / S

    else:
        S = sqrt(1 + m33 - m11 -m22) * 2

        w = (m21 - m12) / S
        x = (m13 + m31) / S
        y = (m23 + m32) / S
        z = 0.25 * S

    return np.array([w, x, y, z])

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
        r = t44as33 + 1
        s = 0.5 / sqrt(r)
        w = s * r
        x = (m23 - m32) * s
        y = (m31 - m13) * s
        z = (m12 - m21) * s
        # w = (m23 - m32) * s
        # x = (m31 - m13) * s
        # y = (m12 - m21) * s
        # z = s * r

    elif m11 > m22 and m11 > m33:
        r = 1.0 - t44as33 + 2 * m11
        #r = m11  - m22 - m33 + 1
        s = 0.5 / sqrt(r)
        w = (m23 - m32) * s
        x = s * r
        y = (m12 + m21) * s
        z = (m13 + m31) * s
        # w = s * r
        # x = (m12 + m21) * s
        # y = (m13 + m31) * s
        # z = (m23 - m32) * s

    elif m22 > m33:
        r = 1.0 - t44as33 + 2 * m22
        #r = 1 - m11 + m22 - m33
        s = 0.5 / sqrt(r)
        w = (m31 - m13) * s
        x = (m12 + m21) * s
        y = s * r
        z = (m23 + m32) * s
        # w = (m12 + m21) * s
        # x = s * r
        # y = (m23 + m32) * s
        # z = (m31 - m13) * s

    else:
        r = 1.0 - t44as33  + 2 * m33
        #r = 1.0 - m11 - m22 + m33
        s = 0.5 / sqrt(r)
        w =  (m12 - m21) * s
        x = (m31 + m13) * s
        y = (m23 + m32) * s
        z = s * r
        # w = (m31 + m13) * s
        # x = (m23 + m32) * s
        # y = s * r
        # z = (m12 - m21) * s

    return np.array([w,x,y,z])

def _quatarc(quat1, quat2, step):
    from math import acos, sqrt, sin

    coshaltheta = quat1[0] * quat2[0] + quat1[1] * quat2[1] + quat1[2] * quat2[2] + quat1[3] * quat2[3]

    if coshaltheta < 0:
        quat2 = quat2 * [-1, -1, -1, 1]
        coshaltheta = -coshaltheta

    if abs(coshaltheta) >= 1:
        quat = quat1
        return quat

    halftheta =  acos(coshaltheta)
    sinhalftheta = sqrt(1- coshaltheta**2)

    if abs(sinhalftheta) < 0.001:
        quat = np.array([quat1[0] * 0.5 + quat2[0] * 0.5,
                         quat1[1] * 0.5 + quat2[1] * 0.5,
                         quat1[2] * 0.5 + quat2[2] * 0.5,
                         quat1[3] * 0.5 + quat2[3] * 0.5])
        return quat

    ratio1 = sin(1-step) * halftheta / sinhalftheta
    ratio2 = sin(step * halftheta) / sinhalftheta

    quat = np.array([quat1[0] * ratio1 + quat1[0] * ratio2,
                     quat1[1] * ratio1 + quat1[1] * ratio2,
                     quat1[2] * ratio1 + quat1[2] * ratio2,
                     quat1[3] * ratio1 + quat1[3] * ratio2])

    return quat

def quatarc(quat1, quat2, step):
    from math import sqrt, acos, sin, fabs, degrees, atan, atan2, cos

    qdot = np.dot(quat1, quat2)

    if qdot > 0.9999:
        return quat2
    elif qdot < 0:
        quat1 = quat1 * -1

    theta = acos( np.dot(quat1, quat2)/sqrt(np.dot(quat1, quat1) * np.dot(quat2, quat2) ) )

    quat = np.add( quat1 * sin(theta * (1-step))/sin(theta), quat2 * sin(theta * step)/sin(theta) )




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

    new_euler = euler

    for i in range(diff_euler.shape[0]):
        if diff_euler[i] > pi:
            new_euler[i] = euler[i]-2*pi
        elif diff_euler[i] < -pi:
            new_euler[i] = euler[i] + 2 * pi
    return np.array(new_euler)

def _invQuaternion(quat):
    quatInv = quat * [1,-1,-1,-1]

    return quatInv

def _normalizeQuat(q):
    from math import sqrt
    d = sqrt(sum( [ i**2 for i in q] ))
    return np.array([i/d for i in q])