from htmd.vmdviewer import VMD, getCurrentViewer
from nglview import NGLWidget,HTMDTrajectory
from tempfile import NamedTemporaryFile
import numpy as np
import os
from glob import glob
import logging
import shutil
import moviepy.editor as mpy


logger = logging.getLogger(__name__)

class Scene:

    def __init__(self, id=None, reps=None, frame=None, delay=1, log=True):

        self.id = id
        if frame == None and log:
            logger.warning('Frame argument not passed. The scene will be applied at the initial frame')
            frame = 0
        self.frame = frame
        self.delay = delay

        self.rollEnabled = False

        self.roll_params = {'selection': 'protein and name CA',
                            'axis':'z',
                            'degree':10,
                            'steps':36,
                            'delay':0.2}

        # vmd
        self.representations = reps
        self.rotate_matrix = None
        self.center_matrix = None
        self.scale_matrix = None
        self.global_matrix = None

        # ngl
        self.orientation_matrix = None
        self.orientation_delta = None

    def setRoll(self, sel='protein and name CA', axis='z', degree=10, steps=36, delay=0.2):
        self.rollEnabled = True
        self.roll_params = {'selection':sel,
                         'axis':axis,
                         'degree':degree,
                         'steps':steps,
                         'delay':delay}


class Transaction:

    def __init__(self, begScene, endScene, scenes=None, id=None, frame=None, delay=1, log=True):

        if frame == None and log:
            logger.warning('Frame argument not passed. The transaction will be applied at the initial frame')
            frame = 0
        self.id = id
        self.frame = frame
        self.delay = delay

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

        self._currentRepresentations = None

        self._initStage()

    def _initStage(self):
        mol = self.mol
        viewer = self.viewer

        viewer.setBackground(self.background)
        viewer.loadMolecule(mol)

    def _getScene(self, frame, id, clonereps, fromscene):
        # repr current, previous
        viewer = self.viewer
        mol = self.mol

        scene = Scene(id=id,  reps=mol.reps, frame=frame)

        r_mat, c_mat, s_mat, g_mat, o_mat = viewer.getMatrices()

        scene.rotate_matrix = r_mat
        scene.center_matrix = c_mat
        scene.scale_matrix = s_mat
        scene.global_matrix = g_mat
        scene.orientation_matrix = o_mat

        if clonereps == True:
            ref_scene = self.scenes[fromscene]
            reps = ref_scene.representations
        else:
            reps = viewer.getReps(mol)

        scene.representations = reps

        return scene

    def saveScene(self, frame=None, clonereps=False, fromscene=-1):
        _id = len(self.scenes)

        scene = self._getScene(frame, _id, clonereps, fromscene)

        self.scenes.append(scene)

    def retrieveScene(self, scene):

        if isinstance(scene, Scene):
            scene = scene

        else:
            if scene >= len(self.scenes):
                raise IndexError('The scene with that sceneId does not exists')
            scene = self.scenes[scene]

        viewer = self.viewer

        change = False
        if self._currentRepresentations is None:
            change = True
            self._currentRepresentations = scene.representations
        else:
            curr_rep = [[r.sel, r.style, r.color, r.material, r.selupdate] for r in
                        self._currentRepresentations.replist]
            new_rep = [[r.sel, r.style, r.color, r.material, r.selupdate] for r in scene.representations.replist]
            if curr_rep != new_rep:
                change = True
                self._currentRepresentations = scene.representations

        viewer.renderScene(scene, updateReps=change)

    def record(self, moviename='mymovie.avi', outdir="imagesVideo", skip=1, overwrite=False, fps=None, spf=None):
        viewer = self.viewer

        extension = os.path.splitext(moviename)[-1]
        if len(extension) == 0:
            moviename = moviename + '.avi'

        elif extension != '.avi':
            raise ValueError('The movie file extension should be \'avi\'.')

        if os.path.isdir(outdir):
            if not overwrite:
                raise IsADirectoryError('The directory exists, change the name or set as True the overwrite argument')
            else:
                shutil.rmtree(outdir)
        os.mkdir(outdir)

        viewer._setUpRecording()

        self.play(skip=skip, _record=True, _outdir=outdir)
        images = glob(os.path.join(outdir, '*.png'))

        # TODO duration, fps, spf
        clip = mpy.ImageSequenceClip(sorted(images), fps=24)#, durations=np.ones(len(images))*0.2)
        clip.write_videofile(moviename, fps=24, codec='png')

    def _playTraj(self, viewer, numFrames, timeline, skip, delay, _record, _outdir):
        from time import sleep
        from tqdm import trange


        play_frames = list(range(numFrames))

        play_frames = [f for i, f in enumerate(play_frames) if i % skip == 0]

        for el in timeline:
            place_idx = play_frames.index(el.frame) + 1
            while not isinstance(play_frames[place_idx], int):
                place_idx += 1
            play_frames.insert(place_idx, el)



        print(play_frames)

        for  i in trange(len(play_frames)):
            el = play_frames[i]
            if _record:
                nf = "%06d" % len(glob(os.path.join(_outdir, '*.png')))
                fname = os.path.join(_outdir, "image.{}".format( nf))

            if isinstance(el, int):
                viewer.goToFrame(el*skip)
                if _record:
                    viewer.render(fname)
                sleep(delay)
            elif isinstance(el, Scene):
                print('rendering scene')
                self.retrieveScene(el)
                if _record:
                    viewer.render(fname)
                #TODO check it
                if el.rollEnabled:
                    viewer.applyRoll(el, _record, _outdir)

                sleep(el.delay)
            elif isinstance(el, Transaction):

                self.playTransaction(el, _record, _outdir)
                sleep(el.delay)

    def _playStructure(self, viewer, timeline):
        pass

    def play(self,  delay=0.3, skip=1, _record=False, _outdir=None):

        viewer = self.viewer

        if isinstance(viewer, NGLWidget):
            print('Not available. At the moment it does not work')
            return

        mol = self.mol
        if len(self.timeline) == 0:
            self.autoTimeline()
        timeline = self.timeline

        numFrames = mol.numFrames

        if numFrames > 1:
            self._playTraj(viewer, numFrames, timeline, skip, delay, _record, _outdir)
        else:
            if len(timeline) == 0:
                raise Exception('The Molecule object is a single frames with not scenes created.')
            self._playStructure(viewer, timeline)


    def playTransaction(self, transaction, _record=False, _outdir=None):

        if _record and _outdir is None:
            raise ValueError('_outdir argument required if _record is True')

        viewer = self.viewer

        for scene in transaction.scenes:
            self.retrieveScene(scene)

            if _record:
                n = "%06d" % len(glob(os.path.join(_outdir, '*.png')))
                fname = os.path.join(_outdir, "image.{}".format(n))
                viewer.render(fname)

    def autoTimeline(self):
        from collections import Counter

        tmp_timeline = self.scenes + self.animations

        tmp_timeline = sorted(tmp_timeline, key= lambda k: k.frame)
        tmp_breakpoints = {e: e.frame for e in tmp_timeline}

        breakpoints_counts = Counter(tmp_breakpoints.values())

        timeline = []
        for k in breakpoints_counts.keys():
            counts = breakpoints_counts[k]
            if counts > 1:
                elements = [el for el, fr in tmp_breakpoints.items() if fr == k]
                scenes = [el for el in elements if isinstance(el, Scene)]
                scenes = sorted(scenes, key= lambda k: k.id)
                timeline.extend(scenes)
                transactions = [el for el in elements if isinstance(el, Transaction)]
                transactions = sorted(transactions, key=lambda k: k.id)
                timeline.extend(transactions)

            else:
                el = [el for el, fr in tmp_breakpoints.items() if fr == k][0]
                timeline.append(el)

        #print("Timeline:", [(el, el.frame, el.id) for el in timeline])
        self.timeline = timeline

    def transaction(self, begSceneId, endSceneId, numsteps=50, frame=None):

        _id = len(self.animations)

        begScene = self.scenes[begSceneId]
        endScene = self.scenes[endSceneId]

        viewer = self.viewer

        listScenes = viewer.setUpTransaction(begScene, endScene, numsteps, frame)

        T = Transaction(begScene, endScene, listScenes, id=_id, frame=frame)
        self.animations.append(T)

    def _getViewer(self, viewer):

        if isinstance(viewer, VMD):
            if  viewer.completed():
                return VMDviewer(getCurrentViewer())
            return VMDviewer(viewer)

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


    elif m11 > m22 and m11 > m33:
        r = 1.0 - t44as33 + 2 * m11
        s = 0.5 / sqrt(r)
        w = (m23 - m32) * s
        x = s * r
        y = (m12 + m21) * s
        z = (m13 + m31) * s

    elif m22 > m33:
        r = 1.0 - t44as33 + 2 * m22
        s = 0.5 / sqrt(r)
        w = (m31 - m13) * s
        x = (m12 + m21) * s
        y = s * r
        z = (m23 + m32) * s

    else:
        r = 1.0 - t44as33  + 2 * m33
        s = 0.5 / sqrt(r)
        w =  (m12 - m21) * s
        x = (m31 + m13) * s
        y = (m23 + m32) * s
        z = s * r

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
    from math import sqrt, acos, sin

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


class VMDviewer:

    def __init__(self, viewer):

        self.viewer = viewer

        self._initStage()

    def _initStage(self):
        cmd_script ="""display rendermode GLSL
display height 3.0"""

        self.runCmd(cmd_script)

    def setBackground(self, bgcolor):
        cmd = 'color Display Background {}'.format(bgcolor)
        self.runCmd(cmd)

    def loadMolecule(self, mol):
        viewer = self.viewer
        viewer.loadMol(mol)

    def runCmd(self, cmd):

        viewer = self.viewer
        viewer.send(cmd)

    def getMatrices(self):
        rotate_matrix = self._matrixFromVMD('rotate')
        center_matrix = self._matrixFromVMD('center')
        scale_matrix = self._matrixFromVMD('scale')
        global_matrix = self._matrixFromVMD('global')

        return rotate_matrix, center_matrix, scale_matrix, global_matrix, None

    def _matrixToVMD(self, scene, matrixtype):

        mat = getattr(scene, "{}_matrix".format(matrixtype))

        mat = str(mat.tolist())
        mat = mat.replace('[', '{').replace(']', '}').replace(',', '')

        return mat

    def _matrixFromVMD(self, matrixtype):
        viewer = self.viewer

        outputFile = NamedTemporaryFile(delete=False).name

        cmd = 'set R [molinfo top get {}_matrix]'.format(matrixtype)

        self.runCmd(cmd)
        self._writeTclOutput( 'R', outputFile)

        f = open(outputFile, 'r')
        txt_Matrix = f.read().replace('{', '').replace('}', '')
        matrix = np.array(txt_Matrix.split(), dtype=float)
        matrix = matrix.reshape(4, 4)

        return matrix


    def getReps(self, mol):
        from htmd.molecule.molecule import Representations

        outputFile = NamedTemporaryFile(delete=False).name

        script_cmd = """set numreps [molinfo top get numreps]
set repslist {}

for {set i 0} {$i < $numreps} {incr i} {
    puts $i
    set replist [molinfo top get "{rep $i} {selection $i} {color $i} {material $i}"]
    set stringreps [join $replist ", " ]
    lappend repslist $stringreps
}
set stringrepslist [ join $repslist "\n"]
"""
        self.runCmd(script_cmd)
        r = Representations(mol)

        self._writeTclOutput('stringrepslist', outputFile)

        f = open(outputFile, 'r')
        list_reps = f.readlines()
        for rep in list_reps:
            l = [_.strip() for _ in rep.strip().split(',')]
            r.add(l[1], l[0], l[2], l[3])
        return r


    def _writeTclOutput(self, outputVariable, outputFile):
        script_cmd ="""set Fname [open {} "w"]
puts -nonewline $Fname ${}
close $Fname
""".format(outputFile, outputVariable)

        self.runCmd(script_cmd)

    def renderScene(self, scene, updateReps=False):

        viewer = self.viewer

        rot_mat = self._matrixToVMD(scene, 'rotate')
        cent_mat = self._matrixToVMD(scene, 'center')
        scale_mat = self._matrixToVMD(scene, 'scale')
        glob_mat = self._matrixToVMD(scene, 'global')

        script_cmd = 'molinfo top set {{rotate_matrix center_matrix scale_matrix global_matrix}} ' \
                    '{{{} {} {} {}}}'.format(rot_mat, cent_mat, scale_mat, glob_mat)

        self.runCmd(script_cmd)

        if updateReps:
            cmd = "set nreps [molinfo top get numreps]\nfor {set i 0} {$i < $nreps } {incr i}  {mol delrep top 0}"
            self.runCmd(cmd)
            scene.representations._repsVMD(viewer)

    def goToFrame(self, frame):

        cmd = "animate goto {}".format(frame)

        self.runCmd(cmd)

    def _setUpRecording(self):
            cmd_script = """global env
set Arch [vmdinfo arch]
set vmdEnv $env(VMDDIR)
set molID [molinfo top]
set tach "$vmdEnv/tachyon_$Arch"

"""
            self.runCmd(cmd_script)

    def applyRoll(self, scene, _record, _outdir):
        from htmd.rotationmatrix import rotationMatrix
        from math import radians
        from time import sleep
        from argparse import Namespace

        _axis = {'x': 0, 'X': 0,
                 'y': 1, 'Y': 1,
                 'z': 2, 'Z': 2}

        _choices = np.unique(list(_axis.keys()) + list(_axis.values()))

        params = Namespace(**scene.roll_params)

        if params.axis not in _choices:
            raise ValueError('The axis should be a valid one: {}'.format(_choices))

        axis = params.axis if isinstance(params.axis, int) else _axis[params.axis]

        rotate_matrix = self._matrixFromVMD('rotate')

        transpose = np.transpose(rotate_matrix)
        rotate_axis = transpose[axis][:3]
        rotate_matrixAxis = rotationMatrix(rotate_axis, radians(params.degree))
        y = np.zeros((3, 1))
        x = np.zeros((1, 4))
        x[-1][-1] = 1
        rotate_matrixAxis = np.append(rotate_matrixAxis, y, axis=1)
        rotate_matrixAxis = np.append(rotate_matrixAxis, x, axis=0)

        cmd_script = """set old_center [molinfo top get center]
    set sel [atomselect top "{}"]
    set new_center [measure center $sel]
    molinfo top set center  [list $new_center]
                """.format(params.selection)

        self.runCmd(cmd_script)

        for i in range(params.steps):
            rotate_matrix = self._matrixFromVMD( 'rotate')

            rotate_matrix_new = np.matmul(rotate_matrixAxis, rotate_matrix)
            mat = str(rotate_matrix_new.tolist())
            mat = mat.replace('[', '{').replace(']', '}').replace(',', '')
            cmd = "molinfo top set rotate_matrix {{ {} }}".format(mat)
            self.runCmd(cmd)
            if _record:
                n = "%06d" % len(glob(os.path.join(_outdir, '*.png')))
                fname = os.path.join(_outdir, "image.{}".format(n))
                self.render(fname)
            sleep(params.delay)

    def setUpTransaction(self, begScene, endScene, numsteps, frame):

        diff_center = endScene.center_matrix - begScene.center_matrix
        diff_scale = endScene.scale_matrix - begScene.scale_matrix
        diff_global = endScene.global_matrix - begScene.global_matrix

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

            sc = Scene(frame=frame)
            sc.rotate_matrix = current_rotateMatrix
            sc.center_matrix = current_centerMatrix
            sc.scale_matrix = current_scaleMatrix
            sc.global_matrix = current_globalMatrix
            sc.representations = begScene.representations

            list_scenes.append(sc)
        return list_scenes

    def render(self, fname):
        # working for vmd only
        viewer = self.viewer

        cmd_script = """set fname {}.tga
render Tachyon $fname "$tach" -aasamples 12 %s -format TARGA -o %s

        """.format(fname)
        viewer.send(cmd_script)
        os.system('convert {}.tga {}.png'.format(fname, fname))
        os.remove('{}.tga'.format(fname))


##### NGL stuffs
# init and load
#  elif initReps and isinstance(viewer, NGLWidget):
        #     traj = HTMDTrajectory(mol)
        #     viewer.add_trajectory(traj)
        #     mol.reps._repsNGL(viewer)

# if isinstance(viewer, NGLWidget):
        #     orientation_matrix = self._matrixFromNGL()
        #     scene.orientation_matrix = orientation_matrix


        # elif isinstance(viewer, NGLWidget):
        #     _class = viewer.control
        #     method = 'orient'
        #     orientation_matrix = self._matrixToNGL(scene)
        #
        #     self._updateView(_class, method, orientation_matrix)
        #     viewer.sync_view()

            #viewer._set_camera_orientation(orientation_matrix)
            #viewer.sync_view()

# if isinstance(viewer, NGLWidget):
#     _class = viewer.control
#     method = 'apply_matrix'
#     mat = scene.orientation_delta
#     args = np.concatenate(mat).tolist()
#     self._updateView(_class, method, args)
#     viewer.sync_view()

# def _transactionNGL(self, begScene, endScene, numsteps, frame):
    #
    #     current_orientation = begScene.orientation_matrix
    #
    #     beg_orientationQuat = _normalizeQuat(matrixToQuaternion(begScene.orientation_matrix))
    #     end_orientationQuat = _normalizeQuat(matrixToQuaternion(endScene.orientation_matrix))
    #
    #     stepsize = 1 / numsteps
    #     current_deltaQuat = quatarc(beg_orientationQuat, end_orientationQuat, stepsize)
    #
    #     current_deltaMatrix = quaternionToMatrix(current_deltaQuat)
    #
    #     list_scenes = []
    #     for i in range(numsteps+1):
    #         sc = Scene(frame=frame)
    #         sc.orientation_delta = current_deltaMatrix
    #
    #         list_scenes.append(sc)
    #
    #     return list_scenes


    # def _matrixToNGL(self, scene):
    #
    #     mat = scene.orientation_matrix
    #     mat = np.concatenate(mat).tolist()
    #
    #     return mat


    # def _matrixFromNGL(self):
    #     viewer = self.viewer
    #
    #     rot_matrix = np.array(viewer._camera_orientation).reshape(4,4)
    #     return rot_matrix

