
# This file is part of the Coriolis Software.
# Copyright (c) Sorbonne Université 2026-2026, All Rights Reserved
#
# +-----------------------------------------------------------------+
# |                   C O R I O L I S                               |
# |      C u m u l u s  -  P y t h o n   T o o l s                  |
# |                                                                 |
# |  Author      :                           Hongda YANG            |
# |  E-mail      :                  laskcl4116@gmail.com            |
# | =============================================================== |
# |  Python      :       "./plugins/block/nldm.py"                  |
# +-----------------------------------------------------------------+

"""
Minimal reader of the NLDM data of a Liberty file: pins capacitances,
max_capacitance and the delay/transition tables of the timing arcs.
Cells are parsed only when requested. Tables are stored with the input
transition as first index and the output load as second, whatever the
order of the variables of their lu_table_template.

The values are in the units of the Liberty file (usually ns and pF).
"""

import re


class NldmCell ( object ):

    def __init__ ( self, name ):
        self.name   = name
        self.pins   = {}      # pin name -> input capacitance.
        self.maxCap = {}      # pin name -> max_capacitance (outputs).
        self.tables = {}      # (relatedPin, pin, kind) -> (index1, index2, values).

    def lookup ( self, kind, inputSlew, load, relatedPin=None, pin=None ):
        """
        Return the value of the ``kind`` table ('cell_rise', 'cell_fall',
        'rise_transition' or 'fall_transition') for the given input transition
        and output load. If the arc is not specified, the first one is used.
        """
        for (rpin, opin, k), table in self.tables.items():
            if k != kind: continue
            if relatedPin is not None and rpin != relatedPin: continue
            if pin        is not None and opin != pin       : continue
            return interpolate( table, inputSlew, load )
        raise KeyError( '{}: no "{}" table'.format(self.name, kind) )


def interpolate ( table, x, y ):
    """Bilinear interpolation, linear extrapolation outside the table."""
    index1, index2, values = table
    def bracket ( index, v ):
        if len(index) < 2: return 0
        for i in range( len(index)-2 ):
            if v < index[i+1]: return i
        return len(index)-2
    i = bracket( index1, x )
    j = bracket( index2, y )
    if len(index1) < 2 and len(index2) < 2: return values[0][0]
    if len(index1) < 2:
        ty = (y - index2[j]) / (index2[j+1] - index2[j])
        return (1-ty)*values[0][j] + ty*values[0][j+1]
    if len(index2) < 2:
        tx = (x - index1[i]) / (index1[i+1] - index1[i])
        return (1-tx)*values[i][0] + tx*values[i+1][0]
    tx = (x - index1[i]) / (index1[i+1] - index1[i])
    ty = (y - index2[j]) / (index2[j+1] - index2[j])
    return (1-tx)*(1-ty)*values[i][j]   + tx*(1-ty)*values[i+1][j] \
         + (1-tx)*ty    *values[i][j+1] + tx*ty    *values[i+1][j+1]


class NldmLibrary ( object ):

    TableKinds = ( 'cell_rise', 'cell_fall', 'rise_transition', 'fall_transition' )

    def __init__ ( self, path ):
        self.path  = path
        with open( path ) as fd:
            self.text = fd.read()
        self.cells = {}
        self.loadFirst = set()
        for name, body in self._subgroups( self.text, 'lu_table_template' ):
            variable1 = re.search( r'variable_1\s*:\s*(\w+)', body )
            if variable1 and 'capacitance' in variable1.group(1):
                self.loadFirst.add( name )

    @staticmethod
    def _group ( text, start ):
        """Return the body of the {...} group opened after ``start``."""
        i     = text.index( '{', start )
        depth = 0
        for j in range( i, len(text) ):
            if   text[j] == '{': depth += 1
            elif text[j] == '}':
                depth -= 1
                if depth == 0: return text[i+1:j]
        raise ValueError( 'NldmLibrary: unbalanced group in "{}".'.format(text[start:start+40]) )

    @staticmethod
    def _subgroups ( text, kind ):
        """Yield (argument, body) of the direct or nested ``kind`` groups of ``text``."""
        for m in re.finditer( r'\b{}\s*\(\s*"?([^")]*)"?\s*\)\s*\{{'.format(kind), text ):
            yield m.group(1).strip(), NldmLibrary._group( text, m.start() )

    @staticmethod
    def _floats ( text ):
        return [ float(v) for v in text.split(',') ]

    def getCell ( self, name ):
        """Return the NldmCell ``name``, or None if the library does not have it."""
        if name in self.cells: return self.cells[name]
        m = re.search( r'\bcell\s*\(\s*"?{}"?\s*\)'.format(re.escape(name)), self.text )
        if not m:
            self.cells[name] = None
            return None
        cell = NldmCell( name )
        body = self._group( self.text, m.start() )
        for pinName, pinBody in self._subgroups( body, 'pin' ):
            cap = re.search( r'(?m)^\s*capacitance\s*:\s*([\d.eE+-]+)', pinBody )
            cell.pins[ pinName ] = float(cap.group(1)) if cap else 0.0
            maxCap = re.search( r'max_capacitance\s*:\s*([\d.eE+-]+)', pinBody )
            if maxCap: cell.maxCap[ pinName ] = float(maxCap.group(1))
            for _, timing in self._subgroups( pinBody, 'timing' ):
                related = re.search( r'related_pin\s*:\s*"?(\w+)"?', timing )
                related = related.group(1) if related else None
                for kind in NldmLibrary.TableKinds:
                    for template, table in self._subgroups( timing, kind ):
                        values = re.search( r'values\s*\((.*?)\)\s*;', table, re.S ).group(1)
                        values = [ self._floats(row) for row in re.findall( r'"([^"]*)"', values ) ]
                        index1 = re.search( r'index_1\s*\(\s*"([^"]*)"', table )
                        index2 = re.search( r'index_2\s*\(\s*"([^"]*)"', table )
                        index1 = self._floats(index1.group(1)) if index1 else [ 0.0 ]
                        index2 = self._floats(index2.group(1)) if index2 else [ 0.0 ]
                        if template in self.loadFirst:
                            index1, index2 = index2, index1
                            values = [ list(column) for column in zip(*values) ]
                        cell.tables[ (related, pinName, kind) ] = (index1, index2, values)
        self.cells[name] = cell
        return cell
