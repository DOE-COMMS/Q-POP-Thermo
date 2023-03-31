# Import Modules
import sympy as sp
import numpy as np
import scipy.optimize as sciop
import itertools as itools
import textwrap
import json, sys
import matplotlib.pyplot as plt

## Class Codes

class terminalPrint:

    def printLine(symbol='*'):
        print(symbol*70)

    def printNamedLine(name, symbol='*'):
        colLength = 70
        nameLength = len(name)
        print('{}{}{}'.format(symbol*int((colLength-nameLength)/2), name, symbol*int((colLength-nameLength)/2)))

    def updateScreen(energyObject):

        printString = ''
        values = postProcess.retrieveParams(energyObject)
        params = energyObject.sweepParams

        for param, varValue in zip(params, values):
            printString += 'Completed: {}: {}; '.format(param, varValue)

        print(printString)

    def printVector(parameters, values):

        terminalPrint.printLine()
        terminalPrint.printNamedLine('Vector')
        for param, varValue in zip(parameters, values):
            print('{}: {}'.format(param, varValue))

        terminalPrint.printLine()

    def printTensor(values, type='Dielectric'):

        terminalPrint.printLine(symbol='-')
        terminalPrint.printNamedLine('{} Tensor'.format(type), symbol='-')
        for iRow in values:
            printString = ''
            for iCol in iRow:
                printString += '{}, '.format(round(iCol,4))
            print(printString)

        terminalPrint.printLine(symbol='-')

    def printIntroMessage():

        terminalPrint.printLine(symbol='-')
        print("     _______ _")
        print("    |__   __| |")
        print(" _   _ | |  | |__   ___ _ __ _ __ ___   ___")
        print("| | | || |  | '_ \ / _ \ '__| '_ ` _ \ / _ \ ")
        print("| |_| || |  | | | |  __/ |  | | | | | | (_) |")
        print("| ____||_|  |_| |_|\___|_|  |_| |_| |_|\___/ ")
        print("| |")
        print("|_|")
        print('')
        try:
        	print('Welcome to \u03BC-Thermo')
        except:
        	print('Welcome to \u03BC-Thermo'.encode('utf-8'))
        print('A Program for Calculating Equilibrium Polarization States')
        print('and properties of Ferroelectric materials.')
        terminalPrint.printLine(symbol='-')
        print('')
        try:
        	print('		\u03BC-Thermo')
        except:
        	print('		\u03BC-Thermo'.encode('utf-8'))
        print('Copyright Notice (c) 2020 Jacob Zorn')
        print('Created at Pennsylvania State University')
        print('Version: 1.0.0')
        print('Last Updated: 10/1/2020')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        print('Help and User Information: Can be found in the accompanying')
        try:
        	print('README file and \u03BC-Thermo manual.')
        except:
        	print('README file and \u03BC-Thermo manual.'.encode('utf-8'))
        print('')
        terminalPrint.printLine(symbol='-')
        print('')
        print('This Software is Licensed via the MIT License.')
        print('Information regarding this license can be found in ')
        print('the provided License.txt software provided with this')
        print('software or at www.gitlab.com/lqc-group/mu_thermo')
        print('')
        terminalPrint.printLine(symbol='-')
        print('')

    def printOutroMessage():

        print('')
        string_to_print = 'The Calculation has completed. Thank you for using mu-Thermo. The data and/or the image files have been properly exported.'
        for text in (textwrap.wrap(string_to_print,width=60)):
        			print(text)
        terminalPrint.printLine(symbol='-')
        print('')
        try:
        	print('		\u03BC-Thermo')
        except:
        	print('		\u03BC-Thermo'.encode('utf-8'))
        print('Copyright Notice (c) 2020 Jacob Zorn')
        print('Created at Pennsylvania State University')
        print('Version: 1.0.0')
        print('Last Updated: 10/1/2020')
        print('')
        terminalPrint.printLine(symbol='-')
        print("     _______ _")
        print("    |__   __| |")
        print(" _   _ | |  | |__   ___ _ __ _ __ ___   ___")
        print("| | | || |  | '_ \ / _ \ '__| '_ ` _ \ / _ \ ")
        print("| |_| || |  | | | |  __/ |  | | | | | | (_) |")
        print("| ____||_|  |_| |_|\___|_|  |_| |_| |_|\___/ ")
        print("| |")
        print("|_|")

class setupRoutines:

    def voigt2ten_array():
        return np.array([[0,5,4], [5,1,3], [4,3,2]])

    def ten2voigt_array():
        return np.array([[0,1,2,1,2,0],[0,1,2,2,0,1]])

    def rank2Transform(mat, ten, toVoigt=True, strain=False):
        if strain:
            m = np.array([1,1,1,2,2,2])
        else:
            m = np.array([1,1,1,1,1,1])

        if toVoigt:
            t2m = setupRoutines.ten2voigt_array()
            for i in range(6):
                mat[i] = ten[t2m[0,i], t2m[1,i]] * m[i]

        if not toVoigt:
            m2t = setupRoutines.voigt2ten_array()
            for i in range(3):
                for j in range(3):
                    ten[i,j] = mat[m2t[i,j]] / m[m2t[i,j]]

        return mat, ten

    def rank3Transform(mat, ten, toVoigt=True):

        if toVoigt:
            t2m = setupRoutines.ten2voigt_array()
            for i in range(3):
                for j in range(6):
                    mat[i,j] = ten[i, t2m[0,j], t2m[1,j]]

        if not toVoigt:
            m2t = setupRoutines.voigt2ten_array()
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        ten[i,j,k] = mat[i,m2t[j,k]]

        return mat, ten

    def rank4Transform(mat, ten, toVoigt=True, compl=False):

        if compl:
            m = np.array([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                              [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                              [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                              [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                              [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                              [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])
        else:
            m = np.ones((6,6))

        if not toVoigt:
            m2t = setupRoutines.voigt2ten_array()
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            ten[i,j,k,l] = mat[m2t[i,j],m2t[k,l]] / m[m2t[i,j], m2t[k,l]]
        else:
            t2m = setupRoutines.ten2voigt_array()
            for i in range(6):
                for j in range(6):
                    mat[i,j] = ten[t2m[0,i],t2m[1,i],t2m[0,j],t2m[1,j]] * m[m2t[i,j], m2t[k,l]]


        return mat, ten

    def rotMatrix(phi, theta, psi, algo=[3,1,3]):
        rotMat = (sp.rot_axis1, sp.rot_axis2, sp.rot_axis3)
        phi, theta, psi = np.deg2rad(phi), np.deg2rad(theta), np.deg2rad(psi)

        transform = (rotMat[algo[2]-1](psi) *
                    rotMat[algo[1]-1](theta) *
                     rotMat[algo[0]-1](phi))

        return transform

    def createRank4_cubic(a11, a12, a44):

        return sp.tensor.MutableDenseNDimArray([[a11,a12,a12,0,0,0],
                                                  [a12,a11,a12,0,0,0],
                                                  [a12,a12,a11,0,0,0],
                                                  [0,0,0,a44,0,0],
                                                  [0,0,0,0,a44,0],
                                                  [0,0,0,0,0,a44]])

    def replaceItems(totalvars, findvars, replacevars):

        for find, replace in zip(findvars, replacevars):
            totalvars[totalvars.index(find)] = replace

        return totalvars

class setupEquations:

    def elasticEnergy(sym, prop):

        energySum = 0
        for i,j,k,l in itools.product(range(3),range(3),range(3),range(3)):
            energySum += 0.5 * prop[i,j,k,l] * sym[i,j] * sym[k,l]

        return energySum

    def couplingEnergy(sym1, sym2, prop):

        energySum = 0
        for i,j,k,l in itools.product(range(3), range(3), range(3), range(3)):
            energySum += prop[i,j,k,l] * sym1[i,j] * sym2[k] * sym2[l]

        return energySum

    def electricEnergy(field, sym, dielectric):

        energySum = 0
        for i in range(3):
            energySum += field[i] * sym[i]

        for i,j in itools.product(range(3), range(3)):
            energySum += 0.5 * field[i] * field[j] * dielectric[i,j]

        return energySum

    def landauEnergy(aTerms, syms):

        a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123 = aTerms

        Landau = 0
        for p in syms:
            Landau += a1 * p**2 + a11 * p**4 + a111 * p**6 + a1111 * p**8
            for p1 in syms:
                if p != p1:
                    Landau += (a12 * p1**2 * p**2)/2
                    Landau += (a112 * p1**2 * p**4)
                    Landau += (a1112 * p1**2 * p**6)
                    Landau += (a1122 * p1**4 * p**4)/2
        Landau += a123 * syms[0]**2 * syms[1]**2 * syms[2]**2
        Landau += a1123 * (syms[0]**4 * (syms[1]**2 * syms[2]**2))
        Landau += a1123 * (syms[1]**4 * (syms[0]**2 * syms[2]**2))
        Landau += a1123 * (syms[2]**4 * (syms[1]**2 * syms[0]**2))

        return Landau

    def thinfilm(bcs, var, energyObject):

        eqs = [
            energyObject.energy.diff(var[0]) + bcs[0],
            energyObject.energy.diff(var[1]) + bcs[1],
            energyObject.energy.diff(var[2]) + bcs[2]
        ]

        solutions = sp.solve(eqs, var)

        for key in energyObject.potential.keys():
            try:
                energyObject.potential[key] = energyObject.potential[key].subs(solutions)
            except:
                continue

        energyObject.energy = energyObject.energy + bcs[0] * var[0] + bcs[1] * var[1] + bcs[2] * var[2]
        energyObject.energy = energyObject.energy.subs(solutions)
        energyObject.problemvars = setupRoutines.replaceItems(energyObject.problemvars, var, bcs)

        return energyObject

    def multilayer(bcs, var, energyObject):

        eqs = [

        ]

    def legredreTransform(energyObject, stress, strain):

        lT = 0
        eqs = []
        for v, vv in zip(stress, strain):
            lT += v * vv
            eqs.append(energyObject.energy.diff(v) + vv)

        solutions = sp.solve(eqs, stress)
        energyObject.energy = energyObject.energy + lT
        energyObject.energy = energyObject.energy.subs(solutions)

        return energyObject

class problem():

    def setProblem(self, inputDict):
        #Set Defaults
        self.system='Polar'
        self.energy = ''
        self.problemvars = ''
        self.bulk = True
        self.film = False
        self.multilayer = False
        self.outputData = []
        self.Stress = [0,0,0,0,0,0]
        self.Strain = [0,0,0,0,0,0]
        self.Fraction = 0.0
        self.electric = False
        self.elastic = False
        self.oxyocta = False
        self.rotate = False
        self.Temp = 298
        self.postprocess = []
        self.outfile = 'test'
        self.clamped = False
        self.solver = 'newton'
        self.solveError = 1e-6
        self.solveIters = 500
        self.parallelWorkers = 1

        #Set according to input file
        self.system = inputDict['Material']
        if inputDict['System'].lower() == 'film':
            self.bulk = False; self.film = True; self.multilayer = False
            self.Misfit = inputDict['Misfit']
        if inputDict['System'].lower() == 'multilayer':
            self.bulk, self.film, self.multilayer = False, False, True
            self.Misfit = inputDict['Misfit']
        if inputDict['System'].lower() == 'bulk':
            self.bulk, self.film, self.multilayer = True, False, False
            try:
                self.Strain = inputDict['Strain']
                self.clamped = True
            except:
                self.Stress = inputDict['Stress']
                self.clamped = False
        if inputDict['System'].lower() == 'multilayer':
            self.bulk, self.film, self.multilayer = False, False, True
            self.Misfit = inputDict['Misfit']
        self.elastic = inputDict['LElastic']
        self.electric = inputDict['LElectric']
        self.rotate = inputDict['LRotate']
        self.RotAngle = inputDict['RotAngle']
        self.oxyocta = inputDict['LRoto']
        self.Temp = inputDict['Temp']
        self.BackgroundDielectric = inputDict['BackgroundDielectric']
        self.EField = inputDict['ElectricField']
        self.sweepList = inputDict['SweepList']
        self.sweepParams = inputDict['SweepParameters']
        self.postprocess = inputDict['Properties']
        self.outfile = inputDict['OutputFile']
        self.solver = (inputDict['Solver']).lower()
        self.solveError = inputDict['SolverError']
        self.solveIters = inputDict['SolverIterations']
        self.parallelWorkers = inputDict['ParallelWorkers']

class calculate:

    def updateArgs(energyObject, parameter, value):
        if parameter.lower() == 'isotropic misfit':
            energyObject.argsDict[sp.symbols('u11')] = value
            energyObject.argsDict[sp.symbols('u22')] = value
        elif parameter.lower() == 'hydrostatic':
            energyObject.argsDict[sp.symbols('s1')] = value
            energyObject.argsDict[sp.symbols('s2')] = value
            energyObject.argsDict[sp.symbols('s3')] = value
        elif parameter.lower() == 'biaxial clamping':
            energyObject.argsDict[sp.symbols('e1')] = value
            energyObject.argsDict[sp.symbols('e2')] = value
        else:
            energyObject.argsDict[sp.symbols(parameter)] = value

    def calculate(energyObject):

        Args = list(energyObject.argsDict.values())

        if energyObject.bulk and not energyObject.clamped:

            def objFunc(x,a,b):
                return energyObject.lambdaproblem(a*x, b*x, x, *Args)

            try:
                xx0 = energyObject.polarVector[2]
            except:
                xx0 = 0.5

            if energyObject.solver == 'newton':
                tetResult = sciop.minimize(objFunc, xx0, args=(0,0), tol=energyObject.solveError, options={'maxiter': energyObject.solveIters})
                ortResult = sciop.minimize(objFunc, xx0, args=(0,1), tol=energyObject.solveError, options={'maxiter': energyObject.solveIters})
                rhoResult = sciop.minimize(objFunc, xx0, args=(1,1), tol=energyObject.solveError, options={'maxiter': energyObject.solveIters})
            elif energyObject.solver == 'differential evolution':
                tetResult = sciop.differential_evolution(objFunc, [(0,1)], args=(0,0), maxiter=energyObject.solveIters, tol=energyObject.solveError, popsize=50)
                ortResult = sciop.differential_evolution(objFunc, [(0,1)], args=(0,1), maxiter=energyObject.solveIters, tol=energyObject.solveError, popsize=50)
                rhoResult = sciop.differential_evolution(objFunc, [(0,1)], args=(1,1), maxiter=energyObject.solveIters, tol=energyObject.solveError, popsize=50)
            elif energyObject.solver == 'simulated annealing':
                tetResult = sciop.dual_annealing(objFunc, [(1e-10,1)], args=(0,0), maxiter=energyObject.solveIters, x0=np.array([xx0]))
                ortResult = sciop.dual_annealing(objFunc, [(1e-10,1)], args=(0,1), maxiter=energyObject.solveIters, x0=np.array([xx0]))
                rhoResult = sciop.dual_annealing(objFunc, [(1e-10,1)], args=(1,1), maxiter=energyObject.solveIters, x0=np.array([xx0]))

            myIDX = [tetResult.fun, ortResult.fun, rhoResult.fun].index(min(tetResult.fun, ortResult.fun, rhoResult.fun))
            if myIDX == 2:
                polar = np.array([rhoResult.x[0], rhoResult.x[0], rhoResult.x[0]])
            elif myIDX == 1:
                polar = np.array([0, ortResult.x[0], ortResult.x[0]])
            else:
                polar = np.array([0,0,tetResult.x[0]])

            result = [tetResult, ortResult, rhoResult][myIDX]

        if energyObject.film or energyObject.clamped:

            def objFunc(x):
                return energyObject.lambdaproblem(*x, *Args)

            try:
                xx0 = energyObject.polarVector
            except:
                xx0 = np.array([0,0,0])

            if energyObject.solver == 'newton':
                result = sciop.minimize(objFunc, xx0, tol=energyObject.solveError, options={'maxiter':energyObject.solveIters})
            elif energyObject.solver == 'differential evolution':
                result = sciop.differential_evolution(objFunc, [(0,1),(0,1),(0,1)], maxiter=energyObject.solveIters, tol=energyObject.solveError)
            elif energyObject.solver == 'simulated annealing':
                result = sciop.dual_annealing(objFunc, [(0,1), (0,1), (0,1)], maxiter=energyObject.solveIters)

            polar = result.x

        return result, polar

    def parameterSweep(energyObject, sweepList, sweepParams):
        sweeper = sweepList[0]; sweepParam = sweepParams[0]
        sweepList = sweepList[1:]; sweepParams = sweepParams[1:]
        for var in np.arange(sweeper[0], sweeper[1]+sweeper[2]/2, sweeper[2]):
            if len(sweepList) > 0:
                calculate.updateArgs(energyObject, sweepParam, var)
                calculate.parameterSweep(energyObject, sweepList, sweepParams)
            else:
                calculate.updateArgs(energyObject, sweepParam, var)
                result, polarVec = calculate.calculate(energyObject)
                energyObject.energy = result.fun
                energyObject.polarVector = np.abs(polarVec)
                postProcessData = postProcess.postProcess(energyObject)
                energyObject.outputData.append(postProcessData)

class postProcess:

    def setupPostProcess(energyObject):

        header = ''
        for param in energyObject.sweepParams:
            header += '{},'.format(param)
        header += 'p1,p2,p3,energy,'

        for prop in energyObject.postprocess:

            if prop.lower() == 'dielectric':
                energyObject.suscept = sp.Matrix([[energyObject.energy.diff(i,j) for i in energyObject.problemvars[:3]] for j in energyObject.problemvars[:3]])
                energyObject.suscept = sp.lambdify(energyObject.problemvars, energyObject.suscept, 'numpy')

                header += 'k11,k22,k33,k23,k13,k12,'

            if prop.lower() == 'entropy':
                energyObject.entropy = -energyObject.energy.diff(sp.symbols('T'))
                energyObject.entropy = sp.lambdify(energyObject.problemvars, energyObject.entropy, 'numpy')

                header += 'S,'

            if prop.lower() == 'piezoelectric':
                stresses = energyObject.problemvars[5:11]

                if energyObject.film:
                    stresses = setupRoutines.replaceItems(stresses, [*sp.symbols('u11 u22 u12')], [*sp.symbols('s1 s2 s6')])

                if energyObject.clamped:
                    stresses = setupRoutines.replaceItems(stresses, stresses, [*sp.symbols('s1 s2 s3 s4 s5 s6')])

                energyObject.strains = sp.Matrix([[-energyObject.energy.diff(i,j) for i in stresses] for j in energyObject.problemvars[:3]])
                energyObject.strains = sp.lambdify(energyObject.problemvars, energyObject.strains, 'numpy')

                header += 'd11,d21,d31,d12,d22,d32,d13,d23,d33,d14,d24,d34,d15,d25,d35,d16,d26,d36,'

        energyObject.header = header

    def retrieveParams(energyObject):

        values = []
        for parameter in energyObject.sweepParams:
            if parameter.lower() == 'isotropic misfit':
                values.append(energyObject.argsDict[sp.symbols('u11')])
            elif parameter.lower() == 'hydrostatic':
                values.append(energyObject.argsDict[sp.symbols('s1')])
            elif parameter.lower() == 'biaxial clamping':
                values.append(energyObject.argsDict[sp.symbols('e1')])
            else:
                values.append(energyObject.argsDict[sp.symbols(parameter)])
        return values

    def postProcess(energyObject):

        params = postProcess.retrieveParams(energyObject)
        terminalPrint.updateScreen(energyObject)

        data = [*params, *energyObject.polarVector, energyObject.energy]

        for prop in energyObject.postprocess:
            if prop.lower() == 'dielectric':
                suscept = energyObject.suscept(*energyObject.polarVector, *list(energyObject.argsDict.values()))
                suscept = 1 + np.linalg.inv(suscept)/8.85e-12
                data = [*data, suscept[0,0], suscept[1,1], suscept[2,2], suscept[1,2], suscept[0,2], suscept[0,1]]

            if prop.lower() == 'entropy':
                data = [*data, energyObject.entropy(*energyObject.polarVector, *list(energyObject.argsDict.values()))]

            if prop.lower() == 'piezoelectric':
                suscept = np.linalg.inv(energyObject.suscept(*energyObject.polarVector, *list(energyObject.argsDict.values())))
                strains = energyObject.strains(*energyObject.polarVector, *list(energyObject.argsDict.values()))
                dij = np.einsum('jn,ij->in',strains,suscept) * 1e12
                dij = [dij[i,j] for j in range(6) for i in range(3)]
                data = [*data,*dij]

        return data

class Input:

    def defaultInput():

        inputDict = {}

        inputDict['Material'] = 'Polar'
        inputDict['System'] = 'Bulk'
        inputDict['Temp'] = 298
        inputDict['LElastic'] = False
        inputDict['LElectric'] = False
        inputDict['LRotate'] = False
        inputDict['LRoto'] = False
        inputDict['RotAngle'] = [0,0,0]
        inputDict['XF'] = 0.00
        inputDict['Stress'] = [0,0,0,0,0,0]
        inputDict['BackgroundDielectric'] = [1,1,1,1,1,1]
        inputDict['ElectricField'] = [0,0,0]
        inputDict['SweepList'] = [[]]
        inputDict['SweepParameters'] = [[]]
        inputDict['Properties'] = []
        inputDict['OutputFile'] = 'Test'

        return inputDict

    def readInputFile(fileName='muInput.in'):

        inputDict = Input.defaultInput()

        with open(fileName,'r') as jsonFile:
            fileDict = json.load(jsonFile)

        for key in fileDict.keys():
            inputDict[key] = fileDict[key]

        return inputDict

    def readPotentialFile(fileName='muPot.in'):

        #Let's make the default potential dictionary
        potDict = {}
        potDict['a1'] = 0 ; potDict['a11'] = 0; potDict['a12'] = 0
        potDict['a111'] = 0; potDict['a112'] = 0; potDict['a123'] = 0
        potDict['a1111'] = 0; potDict['a1112'] = 0; potDict['a1122'] = 0; potDict['a1123'] = 0
        potDict['b1'] = 0; potDict['b11'] = 0; potDict['b12'] = 0
        potDict['S11'] = 0; potDict['S12'] = 0; potDict['S44'] = 0
        potDict['Q11'] = 0; potDict['Q12'] = 0; potDict['Q44'] = 0
        potDict['M11'] = 0; potDict['M12'] = 0; potDict['M44'] = 0
        potDict['l11'] = 0; potDict['l12'] = 0; potDict['l44'] = 0

        with open(fileName,'r') as jsonFile:
            inputPot = json.load(jsonFile)

        for key in inputPot.keys():
            if key in potDict.keys():
                potDict[key] = sp.sympify(inputPot[key])

        return potDict

def main():

    terminalPrint.printIntroMessage()
    myProblem = problem()

    inputDict = Input.readInputFile()
    myProblem.setProblem(inputDict)

    if myProblem.system.lower() == 'polar':

        terminalPrint.printLine(symbol='-')
        terminalPrint.printNamedLine('Ferroelectric Setup Started', symbol=' ')

        terminalPrint.printNamedLine('Starting Landau Setup ...', symbol=' ')
        p1, p2, p3 = sp.symbols('P1 P2 P3')
        a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123 = sp.symbols('a1 a11 a12 a111 a112 a123 a1111 a1112 a1122 a1123')
        polar = sp.Matrix([p1,p2,p3])
        aTerms = [a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123]
        landau = setupEquations.landauEnergy(aTerms, polar)
        totalEnergy = landau
        myProblem.energy = totalEnergy
        probVars = [*polar, sp.symbols('T'), sp.symbols('XF')]
        terminalPrint.printNamedLine('Ending Landau Setup', symbol=' ')

        #Add elastic contribution
        if myProblem.elastic:
            terminalPrint.printNamedLine('Starting Elastic Setup ...', symbol=' ')
            s11, s12, s44 = sp.symbols('S11 S12 S44')
            q11, q12, q44 = sp.symbols('Q11 Q12 Q44')
            s1, s2, s3, s4, s5, s6 = sp.symbols('s1 s2 s3 s4 s5 s6')
            e1, e2, e3, e4, e5, e6 = sp.symbols('e1 e2 e3 e4 e5 e6')
            u11, u22, u12 = sp.symbols('u11 u22 u12')
            stressVector = sp.Matrix([s1, s2, s3, s4, s5, s6])
            strainVector = sp.Matrix([e1, e2, e3, e4, e5, e6])
            stressVector, stressTensor = setupRoutines.rank2Transform(stressVector, sp.Matrix(np.zeros((3,3))), strain=False, toVoigt=False)
            complianceMatrix = setupRoutines.createRank4_cubic(s11, s12, s44)
            complianceMatrix, complianceTensor = setupRoutines.rank4Transform(complianceMatrix, sp.tensor.MutableDenseNDimArray(np.zeros((3,3,3,3))), compl=True, toVoigt=False)
            electroMatrix = setupRoutines.createRank4_cubic(q11, q12, q44)
            electroMatrix, electroTensor = setupRoutines.rank4Transform(electroMatrix, sp.tensor.MutableDenseNDimArray(np.zeros((3,3,3,3))), compl=False, toVoigt=False)
            elastic = setupEquations.elasticEnergy(stressTensor, complianceTensor)
            electroCoupling = setupEquations.couplingEnergy(stressTensor, polar, electroTensor)
            totalEnergy = myProblem.energy - elastic - electroCoupling
            myProblem.energy = totalEnergy
            if myProblem.clamped:
                myProblem = setupEquations.legredreTransform(myProblem, stressVector, strainVector)
                probVars = [*probVars, *strainVector]
            else:
                probVars = [*probVars, *stressVector]

            terminalPrint.printNamedLine('Ending Elastic Setup', symbol=' ')

        #Add Electric contribution
        if myProblem.electric:
            terminalPrint.printNamedLine('Starting Electric Setup ...', symbol=' ')
            E1, E2, E3 = sp.symbols('E1 E2 E3')
            k11, k22, k33, k12, k23, k13 = sp.symbols('k11 k22 k33 k12 k23 k13')
            eField = sp.Matrix([E1, E2, E3])
            dielecVector = sp.Matrix([k11, k22, k33, k23, k13, k12])
            dielecVector, dielecTensor = setupRoutines.rank2Transform(dielecVector, sp.Matrix(np.zeros((3,3))), strain=False, toVoigt=False)
            electric = setupEquations.electricEnergy(eField, polar, dielecTensor)
            totalEnergy = myProblem.energy - electric
            myProblem.energy = totalEnergy
            probVars = [*probVars, *eField, *dielecVector]
            terminalPrint.printNamedLine('Ending Electric Setup', symbol=' ')

        potDict = Input.readPotentialFile()
        myProblem.problemvars = probVars
        myProblem.potential = potDict

        #Check to see if we need to apply thin film boundary conditions
        if myProblem.film:
            terminalPrint.printNamedLine('Starting Thin Film Boundary Conditions ...', symbol=' ')

            #Check to see if we need to rotate the crystal system
            if myProblem.rotate:
                #Create Rotation Matrix
                rotationMatrix = setupRoutines.rotMatrix(*myProblem.RotAngle)
                newOrientation = np.array(rotationMatrix * sp.Matrix([0,0,1]))
                newOrientation = newOrientation / np.max(np.abs(newOrientation))
                newOrientation = 1 / newOrientation
                print(newOrientation)
                rotatedP = rotationMatrix * sp.Matrix([p1,p2,p3])
                rotatedS = rotationMatrix * sp.Matrix([[s1,s6,s5],[s6,s2,s4],[s5,s4,s3]]) * rotationMatrix.transpose()
                p1f, p2f, p3f, s1f, s2f, s3f, s4f, s5f, s6f = sp.symbols('p1^f, p2^f, p3^f, s1^f, s2^f, s3^f, s4^f, s5^f, s6^f')
                fSolutions = sp.solve([rotatedP[0] - p1f, rotatedP[1] - p2f, rotatedP[2] - p3f,
                                      rotatedS[0,0] - s1f, rotatedS[1,1] - s2f, rotatedS[2,2] - s3f,
                                      rotatedS[0,1] - s6f, rotatedS[0,2] - s5f, rotatedS[1,2] - s4f],
                                     [p1,p2,p3,s1,s2,s3,s4,s5,s6])
                myProblem.energy = myProblem.energy.subs(fSolutions)
                for key in myProblem.potential.keys():
                    try:
                        myProblem.potential[key] = myProblem.potential[key].subs(solutions)
                    except:
                        continue
                myProblem.energy = myProblem.energy.subs(dict(zip([p1f, p2f, p3f, s1f, s2f, s3f, s4f, s5f, s6f],
                                                                 [p1, p2, p3, s1, s2, s3, s4, s5, s6])))

            myProblem = setupEquations.thinfilm([u11, u22, u12], [s1, s2, s6], myProblem)
            terminalPrint.printNamedLine('Ending Thin Film Boundary Conditions', symbol=' ')

        myProblem.energy = myProblem.energy.subs(myProblem.potential)
        myProblem.energy = myProblem.energy.subs(myProblem.potential)

        myProblem.lambdaproblem = sp.lambdify(myProblem.problemvars, myProblem.energy, 'numpy')

        if myProblem.film:
            myProblem.argsDict = dict(zip(myProblem.problemvars[3:],
                                         [myProblem.Temp, myProblem.Fraction,
                                         *myProblem.Misfit, *myProblem.EField, *myProblem.BackgroundDielectric]))

        if myProblem.bulk:
            if myProblem.clamped:
                myProblem.argsDict = dict(zip(myProblem.problemvars[3:],
                                         [myProblem.Temp, myProblem.Fraction,
                                         *myProblem.Strain, *myProblem.EField, *myProblem.BackgroundDielectric]))
            else:
                myProblem.argsDict = dict(zip(myProblem.problemvars[3:],
                                         [myProblem.Temp, myProblem.Fraction,
                                         *myProblem.Stress, *myProblem.EField, *myProblem.BackgroundDielectric]))

        terminalPrint.printNamedLine('Starting Post-Processing Setup ...', symbol=' ')
        postProcess.setupPostProcess(myProblem)
        terminalPrint.printNamedLine('Ending Post-Processing Setup ...', symbol=' ')

        terminalPrint.printNamedLine('Ferroelectric Setup Completed', symbol=' ')
        terminalPrint.printLine(symbol='-')

        terminalPrint.printLine(symbol='*')
        terminalPrint.printNamedLine('Beginning Calculations ...', symbol=' ')

        if myProblem.clamped:
            myProblem.bulk, myProblem.film = False, True

        calculate.parameterSweep(myProblem, myProblem.sweepList, myProblem.sweepParams)
        terminalPrint.printNamedLine('End Calculations', symbol=' ')
        terminalPrint.printLine(symbol='*')

        terminalPrint.printLine(symbol='=')
        terminalPrint.printNamedLine('Writing Output to File ...', symbol=' ')
        np.savetxt('{}.out'.format(myProblem.outfile), np.asarray(myProblem.outputData, dtype=object), header=myProblem.header, comments='', delimiter=',')
        terminalPrint.printNamedLine('Output Written to {}.out'.format(myProblem.outfile), symbol=' ')
        terminalPrint.printLine(symbol='=')

    terminalPrint.printOutroMessage()

#Run Code
main()
