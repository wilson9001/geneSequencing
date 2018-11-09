#!/usr/bin/python3
from enum import Enum, auto, IntEnum

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
bandWidth = (2 * MAXINDELS) + 1

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
		self.MaxCharactersToAlign = align_length+1
		results = []
		previous = 'prevhttp://mind.cs.byu.edu/courses/312/projects/project4_files/banded.pngious'
		cost = 'cost'
		topStringAlignment = []
		sideStringAlignment = []

		maxArray = [[{previous: None, cost: math.inf} for x in range(self.MaxCharactersToAlign)] for y in range(self.MaxCharactersToAlign)]
		maxArray1 = [[{previous: None, cost: math.inf} for x in range(bandWidth)] for y in range(self.MaxCharactersToAlign)]

		for row in range(len(sequences)):
			rowResults = []
			for column in range(len(sequences)):

				if(column < row):
					s = {}
				elif column == row:
					s = {'align_cost': -3 * (len(sequences[row]) if len(sequences[row]) < self.MaxCharactersToAlign else self.MaxCharactersToAlign-1), 'seqi_first100': sequences[row][:100], 'seqj_first100': sequences[row][:100]}
					table.item(row, column).setText('{}'.format(s['align_cost']))
					table.repaint()
				else:
###################################################################################################
					# extract sequences for this comparison
					topString = sequences[row]
					sideString = sequences[column]

					topStringAlignmentFinal = None
					sideStringAlignmentFinal = None
					score = None

					if not self.banded:
						# set up matrix
						topStringDashedLen = len(topString)+1 if len(topString)+1 <= self.MaxCharactersToAlign else self.MaxCharactersToAlign
						sideStringDashedLen = len(sideString)+1 if len(sideString)+1 <= self.MaxCharactersToAlign else self.MaxCharactersToAlign

						# set up strings
						topStringDashed = "-" + topString[:topStringDashedLen-1]
						sideStringDashed = "-" + sideString[:sideStringDashedLen-1]

						costMatrix = None

						if topStringDashedLen == self.MaxCharactersToAlign and sideStringDashedLen == self.MaxCharactersToAlign:
							costMatrix = maxArray
						else:
							costMatrix = [[{previous: None, cost: math.inf} for x in range(topStringDashedLen)] for y in range(sideStringDashedLen)]

						# initialize matrix
						costMatrix[0][0][cost] = 0

						for innerColumn in range(1, topStringDashedLen):
							costMatrix[0][innerColumn][previous] = Side.l
							costMatrix[0][innerColumn][cost] = costMatrix[0][innerColumn-1][cost] + INDEL

						for innerRow in range(1, sideStringDashedLen):
							costMatrix[innerRow][0][previous] = Side.t
							costMatrix[innerRow][0][cost] = costMatrix[innerRow-1][0][cost] + INDEL

						for innerRow in range(1, sideStringDashedLen):
							for innerColumn in range(1, topStringDashedLen):

								bestOption = {previous: None, cost: None}

								left = costMatrix[innerRow][innerColumn-1][cost] + INDEL

								diagonal = None

								if topStringDashed[innerColumn] == sideStringDashed[innerRow]:
									diagonal = costMatrix[innerRow-1][innerColumn-1][cost] + MATCH
								else:
									diagonal = costMatrix[innerRow - 1][innerColumn - 1][cost] + SUB

								top = costMatrix[innerRow-1][innerColumn][cost] + INDEL

								if left < diagonal:
									bestOption[previous] = Side.l
									bestOption[cost] = left

								else:
									bestOption[previous] = Side.d
									bestOption[cost] = diagonal

								if top < diagonal:
									bestOption[previous] = Side.t
									bestOption[cost] = top

								costMatrix[innerRow][innerColumn][previous] = bestOption[previous]
								costMatrix[innerRow][innerColumn][cost] = bestOption[cost]

						sideStringIndex = sideStringDashedLen-1
						topStringIndex = topStringDashedLen-1
						currentEntry = costMatrix[sideStringIndex][topStringIndex]

						score = currentEntry[cost]

						topStringAlignment.clear()
						sideStringAlignment.clear()

						while currentEntry[previous] is not None:
							if currentEntry[previous] == Side.l:
								sideStringAlignment.append("-")
								topStringAlignment.append(topStringDashed[topStringIndex])
								topStringIndex -= 1
							elif currentEntry[previous] == Side.d:
								sideStringAlignment.append(sideStringDashed[sideStringIndex])
								topStringAlignment.append(topStringDashed[topStringIndex])
								sideStringIndex -= 1
								topStringIndex -= 1
							else:
								sideStringAlignment.append(sideStringDashed[sideStringIndex])
								topStringAlignment.append("-")
								sideStringIndex -= 1

							currentEntry = costMatrix[sideStringIndex][topStringIndex]

						topStringAlignment.reverse()
						sideStringAlignment.reverse()

						topStringAlignmentFinal = ''.join(topStringAlignment)
						sideStringAlignmentFinal = ''.join(sideStringAlignment)

						# TODO: Create section for amortized later. Don't forget to set score again in new section
					else:
						topStringDashedLen = len(topString) + 1 if len(topString) + 1 <= self.MaxCharactersToAlign else self.MaxCharactersToAlign
						sideStringDashedLen = len(sideString) + 1 if len(sideString) + 1 <= self.MaxCharactersToAlign else self.MaxCharactersToAlign

						diff = topStringDashedLen - sideStringDashedLen

						if abs(diff) > MAXINDELS:
							score = math.inf
							topStringAlignmentFinal = "No Alignment Possible"
							sideStringAlignmentFinal = "No Alignment Possible"
						else:
							# set up strings
							topStringDashed = "-" + topString[:topStringDashedLen - 1]
							sideStringDashed = "-" + sideString[:sideStringDashedLen - 1]

							costMatrix = None

							# create matrix
							if topStringDashedLen == self.MaxCharactersToAlign and sideStringDashedLen == self.MaxCharactersToAlign:
								costMatrix = maxArray1
							else:
								costMatrix = [[{previous: None, cost: math.inf} for x in range(bandWidth)] for y in range(sideStringDashedLen)]

							# fill in top-left corner
							nills = MAXINDELS
							innerRow = 0

							while nills > 0:
								for innerColumn in range(nills):
									costMatrix[innerRow][innerColumn] = None
								innerRow += 1
								nills -= 1

							# fill in bottom-right corner:
							innerRow = sideStringDashedLen-1
							nills = MAXINDELS-diff

							while nills > 0:
								for innerColumn in range(bandWidth-1, (bandWidth-1)-nills, -1):
									costMatrix[innerRow][innerColumn] = None
								innerRow -= 1
								nills -= 1

							startingColumn = MAXINDELS
							innerRow = 0

							# initialize matrix
							costMatrix[0][startingColumn][cost] = 0
							for innerColumn in range(startingColumn+1, bandWidth):
								costMatrix[innerRow][innerColumn][cost] = costMatrix[innerRow][innerColumn-1][cost] + INDEL
								costMatrix[innerRow][innerColumn][previous] = Side.l

							startingColumn -= 1
							topStringMinIndex = 1-MAXINDELS

							for innerRow in range(1, sideStringDashedLen):
								for innerColumn in range(startingColumn, bandWidth):

									if costMatrix[innerRow][innerColumn] is not None:
										bestOption = {previous: None, cost: None}

										left = None

										if innerColumn-1 < 0 or costMatrix[innerRow][innerColumn-1] is None:
											left = math.inf
										else:
											left = costMatrix[innerRow][innerColumn - 1][cost] + INDEL

										diagonal = None
										topStringIndex = topStringMinIndex + startingColumn
										topStringIndex += (innerColumn - startingColumn)

										if topStringIndex > (topStringDashedLen-1):
											topStringIndex = topStringDashedLen-1

										if costMatrix[innerRow - 1][innerColumn] is None:
											diagonal = math.inf
										elif topStringDashed[topStringIndex] == sideStringDashed[innerRow]:
											diagonal = costMatrix[innerRow - 1][innerColumn][cost] + MATCH
										else:
											diagonal = costMatrix[innerRow - 1][innerColumn][cost] + SUB

										top = None

										if innerColumn+1 >= bandWidth:
											top = math.inf
										else:
											top = costMatrix[innerRow - 1][innerColumn + 1][cost] + INDEL

										# determine best option
										if left < diagonal:
											bestOption[previous] = Side.l
											bestOption[cost] = left

										else:
											bestOption[previous] = Side.d
											bestOption[cost] = diagonal

										if top < diagonal:
											bestOption[previous] = Side.t
											bestOption[cost] = top

										costMatrix[innerRow][innerColumn][previous] = bestOption[previous]
										costMatrix[innerRow][innerColumn][cost] = bestOption[cost]

								topStringMinIndex += 1
								startingColumn = startingColumn - 1 if startingColumn > 0 else 0

							sideStringIndex = sideStringDashedLen - 1
							topStringIndex = topStringDashedLen - 1
							innerColumn = bandWidth-1
							innerRow = sideStringIndex

							currentEntry = costMatrix[innerRow][innerColumn]

							while currentEntry is None:
								innerColumn -= 1
								currentEntry = costMatrix[innerRow][innerColumn]

							score = currentEntry[cost]

							topStringAlignment.clear()
							sideStringAlignment.clear()

							while currentEntry[previous] is not None:
								if currentEntry[previous] == Side.l:
									sideStringAlignment.append("-")
									topStringAlignment.append(topStringDashed[topStringIndex])
									topStringIndex -= 1
									innerColumn -= 1
								elif currentEntry[previous] == Side.d:
									sideStringAlignment.append(sideStringDashed[sideStringIndex])
									topStringAlignment.append(topStringDashed[topStringIndex])
									sideStringIndex -= 1
									topStringIndex -= 1
									innerRow -= 1
								else:
									sideStringAlignment.append(sideStringDashed[sideStringIndex])
									topStringAlignment.append("-")
									sideStringIndex -= 1
									innerRow -= 1
									innerColumn += 1

								currentEntry = costMatrix[innerRow][innerColumn]

							topStringAlignment.reverse()
							sideStringAlignment.reverse()

							topStringAlignmentFinal = ''.join(topStringAlignment)
							sideStringAlignmentFinal = ''.join(sideStringAlignment)

					alignment1 = '{}  DEBUG:(seq{}, {} chars,align_len={}{})'.format(topStringAlignmentFinal, row+1, len(sequences[row]), align_length, ',BANDED' if banded else '')
					alignment2 = '{}  DEBUG:(seq{}, {} chars,align_len={}{})'.format(sideStringAlignmentFinal, column+1, len(sequences[column]), align_length, ',BANDED' if banded else '')

###################################################################################################
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(row,column).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()
				rowResults.append(s)
			results.append(rowResults)
		return results


