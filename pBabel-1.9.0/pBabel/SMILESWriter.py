#-------------------------------------------------------------------------------
# . File      : SMILESWriter.py
# . Program   : pDynamo-1.9.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin J. Field (2007-2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""Classes and functions for writing SMILES strings."""

import string

from pCore   import logFile
from pMolecule import PeriodicTable

from SMILESUtilities import CHIRALITYDEFAULTCLASSES, ELEMENTSORGANIC, SMILESConnectivity, VALENCIESORGANIC

#===============================================================================
# . Parameters.
#===============================================================================
# . Bond tokens.
_DOUBLEBONDTOKEN = "="
_NULLBONDTOKEN   = "."
_TRIPLEBONDTOKEN = "#"

#===============================================================================
# . SMILES writer class.
#===============================================================================
class SMILESWriter ( object ):
    """SMILESWriter is the class for writing SMILES strings."""

    def __init__ ( self, **keywordArguments ):
        """Constructor."""
        self.Reset ( )

    def AtomToken ( self, atom ):
        """Return a SMILES token for an atom."""
        # . Check for a reduced atom representation.
        # . Basic checks.
        QREDUCED = ( atom.atomicNumber in ELEMENTSORGANIC ) and ( atom.formalCharge == 0 ) and ( atom.isotope <= 0 )
        # . Check for the expected number of implicit hydrogens (also ensures the atom has standard valence).
        if QREDUCED:
            valence = atom.valence - atom.implicithydrogens
            for v in VALENCIESORGANIC[atom.atomicNumber]:
                hcount = v - valence
                if hcount >= 0: break
            QREDUCED = ( atom.implicithydrogens == hcount )
        # . Check for a reduced chirality representation.
        if atom.chiralityclass is not None:
            QCHIRALITYREDUCED = ( CHIRALITYDEFAULTCLASSES.get ( atom.connections, None ) == atom.chiralityclass ) and ( atom.chiralitynumber <= 4 )
            QREDUCED          = QREDUCED and QCHIRALITYREDUCED
            # . Generate the string.
            if QCHIRALITYREDUCED: ctoken = atom.chiralitynumber * "@"
            else:                 ctoken = "@" + atom.chiralityclass + "{:d}".format ( atom.chiralitynumber )
        else: ctoken = ""
        # . Encode the string.
        if atom.isAromatic: atoken = PeriodicTable.Symbol ( atom.atomicNumber ).lower ( )
        else:              atoken = PeriodicTable.Symbol ( atom.atomicNumber )
        tokens = [ atoken, ctoken ]
        if not QREDUCED:
            tokens[0:0] = "["
            if atom.isotope != 0: tokens[1:1] = "{:d}".format ( isotope )
            if   atom.implicithydrogens == 1: tokens.append ( "H" )
            elif atom.implicithydrogens  > 1: tokens.append ( "H" + "{:d}".format ( atom.implicithydrogens ) )
            if   atom.formalCharge < 0: tokens.append (       "{:d}".format ( atom.formalCharge ) )
            elif atom.formalCharge > 0: tokens.append ( "+" + "{:d}".format ( atom.formalCharge ) )
            tokens.append ( "]" )
        return "".join ( tokens )

    def CreateBranch ( self, current, previous ):
        """Return tokens for a branch."""
        self.QATOMS[current] = True
        tokens   = [ self.AtomToken ( current ) ]
        branches = []
        for bond in current.bonds:
            other = bond.OtherAtom ( current )
            if self.QATOMS.get ( other, False ):
                if other is not previous:
                    if   bond.type.bondOrder == 2: tokens.append ( _DOUBLEBONDTOKEN )
                    elif bond.type.bondOrder == 3: tokens.append ( _TRIPLEBONDTOKEN )
                    crosslink = self.QBONDS.pop ( bond, self.crosslinks + 1 )
                    if crosslink > self.crosslinks:
                        self.crosslinks += 1
                        self.QBONDS[bond] = self.crosslinks
                    tokens.append ( "%" + "{:d}".format ( crosslink ) )
            else:
                branch = []
                if   bond.type.bondOrder == 2: branch.append ( _DOUBLEBONDTOKEN )
                elif bond.type.bondOrder == 3: branch.append ( _TRIPLEBONDTOKEN )
                branch.extend ( self.CreateBranch ( other, current ) )
                branches.append ( branch )
        if len ( branches ) > 0:
            for branch in branches[0:-1]: tokens.extend ( [ "(" ] + branch + [ ")" ] )
            tokens.extend ( branches[-1] )
        return tokens

    def FromSMILESConnectivity ( self, connectivity ):
        """Generate a SMILES from a connectivity."""
        # . Initialization.
        self.Reset ( )
        tokens = []
        # . Loop over isolates.
        for current in connectivity.atoms:
            if not self.QATOMS.get ( current, False ):
                isolate = self.CreateBranch ( current, None )
                if len ( tokens ) > 0: tokens.append ( _NULLBONDTOKEN )
                tokens.extend ( isolate )
        self.Reset ( )
        return "".join ( tokens )

    def Reset ( self ):
        """Reset all attributes."""
        self.QATOMS     = {}
        self.QBONDS     = {}
        self.crosslinks = 0

#===============================================================================
# . Helper functions.
#===============================================================================
def SMILES_FromSystem ( system, **keywordArguments ):
    """Helper function that generates a SMILES from a system."""
    writer       = SMILESWriter ( **keywordArguments )
    connectivity = SMILESConnectivity.FromSystem ( system )
    log = keywordArguments.get ( "log", logFile )
    connectivity.Summary ( log = log )
    smiles       = writer.FromSMILESConnectivity ( connectivity )
    return smiles

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
