from htmd.vmdviewer import VMD, getCurrentViewer
from nglview import NGLWidget,HTMDTrajectory
from tempfile import NamedTemporaryFile
import numpy as np

class Scene:

    def __init__(self, reprs=None):

        # vmd
        self.representations = reprs
        self.rotate_matrix = None
        self.center_matrix = None
        self.scale_matrix = None
        self.global_matrix = None

        # ngl
        self.orientation_matrix = None

class Transaction:

    def __init__(self, begScene, endScene):
        self.begScene = begScene
        self.endScene = endScene

        self.scenes = []

    def appendScene(self, scene):

        self.scenes.append(scene)

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
        elif initReps and isinstance(viewer, NGLWidget):
            traj = HTMDTrajectory(mol)
            viewer.add_trajectory(traj)
            mol.reps._repsNGL(viewer)

    def _getScene(self, repr='current'):
        # repr current, previous
        viewer = self.viewer
        mol = self.mol

        scene = Scene(reprs=mol.reps)

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
        scene = self._getScene()

        self.scenes.append(scene)

    def retrieveScene(self, sceneid):
        try:
            scene = self.scenes[sceneid]
        except:
            raise IndexError('The scene with that sceneId does not exists')

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


    def createAnimation(self, startSceneId, endSceneId):

        start_scene = self.scenes[startSceneId]
        start_rotMat = start_scene.rotate_matrix

        end_scene = self.scenes[endSceneId]
        end_rotMat = end_scene.rotate_matrix

    def transaction(self, begSceneId, endSceneId, numsteps=50):
        self.retrieveScene(begSceneId)

        begScene = self.scenes[begSceneId]
        endScene = self.scenes[endSceneId]

        viewer = self.viewer

        if isinstance(viewer, VMD):
            listScenes = self._transactionVMD(begScene, endScene, numsteps)
        elif isinstance(viewer, NGLWidget):
            pass

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

        for i in range(1, numsteps+1):
            stepratio = stepsize * i
            qarc = quatarc(beg_rotateQuat, end_rotateQuat, stepratio)
            

    def test(self):
        beg_scene = self.scenes[startSceneId]
        end_scene = self.scenes[endSceneId]

        beg_eul = matrixToEuler(beg_scene.rotate_matrix)
        end_eul = matrixToEuler(end_scene.rotate_matrix)

        diff_eul = end_eul - beg_eul

        end_eul = _bestEuler(diff_eul, end_eul)

        diff_eul = end_eul - beg_eul

        diff_center = end_scene.center_matrix - beg_scene.center_matrix
        diff_scale = end_scene.scale_matrix - beg_scene.scale_matrix
        diff_global = end_scene.global_matrix - beg_scene.global_matrix

        _diffScene =  Scene()
        _diffScene.rotate_matrix = diff_eul
        _diffScene.center_matrix = diff_center
        _diffScene.scale_matrix = diff_scale
        _diffScene.global_matrix = diff_global

        beg_rotmat_quat = matrixToQuaternion(beg_scene.rotate_matrix)
        end_rotmat_quat = matrixToQuaternion(end_scene.rotate_matrix)

        beg_centmat_quat = matrixToQuaternion(beg_scene.center_matrix)
        end_centmat_quat = matrixToQuaternion(end_scene.center_matrix)

        beg_scalemat_quat = matrixToQuaternion(beg_scene.scale_matrix)
        end_scalemat_quat = matrixToQuaternion(end_scene.scale_matrix)

        beg_globalmat_quat = matrixToQuaternion(beg_scene.global_matrix)
        end_globalmat_quat = matrixToQuaternion(end_scene.global_matrix)

        self.saveScene()
        now_scene = self.scenes[-1]

        now_centmat = now_scene.center_matrix
        now_scalemat = now_scene.scale_matrix
        now_globmat = now_scene.global_matrix
        import time
        for i in range(50):
            #print(i)
            step = (i+1)/50
            stepsize =  1 /50
            qarc = quatarc(beg_rotmat_quat, end_rotmat_quat, step)
            now_rotmat = quaternionToMatrix(qarc)
            #now_centmat = quaternionToMatrix(quatarc(beg_centmat_quat, end_centmat_quat, step))
            #now_scalemat = quaternionToMatrix(quatarc(beg_scalemat_quat, end_scalemat_quat, step))
            #now_globmat = quaternionToMatrix(quatarc(beg_globalmat_quat, end_globalmat_quat, step))
            now_centmat = np.add(now_centmat, _diffScene.center_matrix * stepsize)
            now_scalemat = np.add(now_scalemat, _diffScene.scale_matrix * stepsize)
            now_globmat = np.add(now_globmat, _diffScene.global_matrix * stepsize)
            rot_mat = str(now_rotmat.tolist())
            rot_mat = rot_mat.replace('[', '{').replace(']', '}').replace(',', '')
            center_mat = str(now_centmat.tolist())
            center_mat = center_mat.replace('[', '{').replace(']', '}').replace(',', '')
            scale_mat = str(now_scalemat.tolist())
            scale_mat = scale_mat.replace('[', '{').replace(']', '}').replace(',', '')
            glob_mat = str(now_globmat.tolist())
            glob_mat = glob_mat.replace('[', '{').replace(']', '}').replace(',', '')
            #print(rot_mat)
            #self.viewer.send('molinfo top set rotate_matrix {{{}}}'.format(rot_mat))
            self.viewer.send('molinfo top set {{rotate_matrix center_matrix scale_matrix global_matrix}} '
                           '{{{} {} {} {}}}'.format(rot_mat, center_mat, scale_mat, glob_mat))
           # time.sleep(0.2)

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
#    print(quat1, quat2)

    # qdot = np.dot(quat1, quat2)
    # print(qdot)
    #
    # if qdot > 0.9999:
    #     return quat2
    # elif qdot < 0:
    #     quat1 = quat1 * -1

    #theta = acos(np.dot(quat1, quat2)/sqrt(np.dot(quat1, quat1) * np.dot(quat2, quat2)))
    #theta = acos(np.dot(quat1, quat2))
    #quat2Conj = _invQuaternion(quat2)
   # print("shape: ", quat1.shape, quat2Conj.shape)
    #print(">>>>> ", np.cross(quat1, quat2Conj))
    # {\displaystyle
    # {} + (a_{1}d_{2}+b_{1}c_{2}-c_{1}b_{2}+d_{1}a_{2})
    # k.} {} + (a_{1}d_{2}+b_{1}c_{2}-c_{1}b_{2}+d_{1}a_{2})
    # k.

    #theta = atan2((quat1 * quat2Conj)[0])
    #print("a: ", degrees(theta), step)
#     cosHalfTheta = np.dot(quat1, quat2)
#     if cosHalfTheta >= 1:
#         return quat1
#
#     halfTheta =  acos(cosHalfTheta)
#     sinHalfTheta = sqrt(1 - cosHalfTheta**2)
#
#     if fabs(sinHalfTheta) < 0.001:
#         quat = np.zeros(4)
#         quat[0] = quat1[0] * 0.5 + quat2[0] * 0.5
#         quat[1] = quat1[1] * 0.5 + quat2[1] * 0.5
#         quat[2] = quat1[2] * 0.5 + quat2[2] * 0.5
#         quat[3] = quat1[3] * 0.5 + quat2[3] ^ 0.5
#         return quat
# #
# # #    quat = np.add( sin(theta * (1-step))/sin(theta) * quat1, sin(theta * step)/sin(theta) * quat2 )
# #     quat = np.add(sin(theta * (1 - step)) * quat1, sin(theta * step) * quat2)  / sin(theta)
# #     print(quat)
# #     print()
#
#     quat = np.zeros(4)
#     ratioA = sin(1-step) * halfTheta / sinHalfTheta
#     ratioB = sin(step * halfTheta) / sinHalfTheta
#
#     quat[0] = quat1[0] * ratioA + quat2[0] * ratioB
#     quat[1] = quat1[1] * ratioA + quat2[1] * ratioB
#     quat[2] = quat1[2] * ratioA + quat2[2] * ratioB
#     quat[3] = quat1[3] * ratioA + quat2[3] * ratioB
#
#     dot = np.dot(quat1, quat2)
#
#     print(dot)
#     if dot < 0:
#         quat1 = -quat1
#         dot = -dot
#
#     threshold = 0.9995
#
#     if dot > threshold:
#         quat = quat1 + step * (quat2 - quat1)
#
#         return quat
#
#     if dot < -1:
#         dot = -1
#     elif dot > 1:
#         dot = 1
#     theta_0 = acos(dot)
#     theta = theta_0 * step
#
#     s0 = cos(theta) - dot * sin(theta)/sin(theta_0)
#     s1 = sin(theta) / sin(theta_0)
#
#     quat = s0 * quat1 + s1 * quat1

    qdot = np.dot(quat1, quat2)

    if qdot > 0.9999:
        return quat2
    elif qdot < 0:
        quat1 = quat1 * -1

    theta = acos( np.dot(quat1, quat2)/sqrt(np.dot(quat1, quat1) * np.dot(quat2, quat2) ) )
    #print(degrees(theta))

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