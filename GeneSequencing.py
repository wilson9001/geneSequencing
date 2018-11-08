#!/usr/bin/python3
from enum import Enum, auto

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class Side(Enum):
	l = auto()
	d = auto()
	t = auto()

class GeneSequencing:

	def __init__( self ):
		pass

	
# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []

		for row in range(len(sequences)):
			rowResults = []
			for column in range(len(sequences)):

				if(column < row):
					s = {}
				else:
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
					# extract sequences for this comparison
					topString = sequences[row]
					sideString = sequences[column]

					# prepend dashes to each string
					topStringDashed = "-" + topString
					sideStringDashed = "-" + sideString

					# set up matrix
					topStringDashedLen = len(topStringDashed)
					sideStringDashedLen = len(sideStringDashed)

					costMatrix = [[math.inf for x in range(topStringDashedLen)] for y in range(sideStringDashedLen)]

					# initialize matrix
					costMatrix[0][0] = {'previous': None, 'cost': 0}

					for innerColumn in range(1, topStringDashedLen):
						costMatrix[0][innerColumn] = {'previous': Side.l, 'cost': costMatrix[0][innerColumn-1] + INDEL}

					for innerRow in range(1, sideStringDashedLen):
						costMatrix[innerRow][0] = {'previous': Side.t, 'cost': costMatrix[innerRow-1][0] + INDEL}

					for innerRow in range(1, sideStringDashedLen):
						for innerColumn in range(1, topStringDashedLen):

							bestOption = {}

							left = costMatrix[innerRow][innerColumn-1] + INDEL
							diagonal = costMatrix[innerRow - 1][innerColumn - 1] + SUB

							if topString[innerColumn] == sideString[row]:
								diagonal = costMatrix[innerRow-1][innerColumn-1] + MATCH

							top = costMatrix[innerRow-1][innerColumn] + INDEL

							if left < diagonal:
								bestOption.previous = Side.l
								bestOption.cost = left
							else:
								bestOption.previous = Side.d
								bestOption.cost = diagonal

							if top < diagonal:
								bestOption.previous = Side.t
								bestOption.cost = top

							costMatrix[innerRow][innerColumn] = dict(bestOption)

					currentEntryRow = sideStringDashedLen-1
					currentEntryColumn = topStringDashedLen-1
					currentEntry = costMatrix[currentEntryRow][currentEntryColumn]

					score = currentEntry.cost

					topStringAlignment = []
					sideStringAlignment = []

					while currentEntry.previous is not None:
						if currentEntry.previous == Side.l:
							sideStringAlignment.append("-")
							topStringAlignment.append(topStringDashed[currentEntryColumn])
							currentEntryColumn -= 1
						elif currentEntry.previous == Side.d:
							sideStringAlignment.append(sideStringDashed[currentEntryRow])
							topStringAlignment.append(topStringDashed[currentEntryColumn])
							currentEntryRow -= 1
							currentEntryColumn -= 1
						else:
							sideStringAlignment.append(sideStringDashed[currentEntryRow])
							topStringAlignment.append("-")
							currentEntryRow -= 1

						currentEntry = costMatrix[currentEntryRow][currentEntryColumn]

					topStringAlignment.reverse()
					sideStringAlignment.reverse()

					# TODO: Create section for amortized later. Don't forget to set score again in new section

					alignment1 = '{}  DEBUG:(seq{}, {} chars,align_len={}{})'.format(''.join(topStringAlignment), row+1, len(sequences[row]), align_length, ',BANDED' if banded else '')
					alignment2 = '{}  DEBUG:(seq{}, {} chars,align_len={}{})'.format(''.join(sideStringAlignment), column+1, len(sequences[column]), align_length, ',BANDED' if banded else '')

					# score = row+column
					# alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(row+1, len(sequences[row]), align_length, ',BANDED' if banded else '')
					# alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(column+1, len(sequences[column]), align_length, ',BANDED' if banded else '')
###################################################################################################					
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(row,column).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()
				rowResults.append(s)
			results.append(rowResults)
		return results


