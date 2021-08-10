""" Purpose of the script is to provide useful functions for converting different
	types of data into different formats.
"""

import numpy, pandas

def listsToFrame(lists, filler='', index=None, columns=None):
	""" Takes of a list of lists, which may or may not be of the same lengths, \
			and adds filler so that they are the same length. Then converts to \
			dataframe with each column containing the inner list information.

	Args:
		lists (list<list<str>>): List of of lists with strings.

		filler (str): Filler value for extra elements to make the lists the \
						same size, always appended to the end of the lists \
						when making the maximum length.

		index (list-like): Values to have on the resulting dataframes rows.

		column (list-like): Value to have on the resulting dataframes columns.

	Returns:
		pandas.DataFrame: Each column contains the list contents; if lists \
							were of unequal length then will have filler.
	"""
	same_length_lists = []
	max_len = max([len(list_set) for list_set in lists])
	for list_set in lists:
		set_len = len(list_set)
		if set_len != max_len:
			list_set = list(list_set) + [filler] * (max_len - set_len)

		same_length_lists.append(list_set)

	lists_frame = pandas.DataFrame(same_length_lists,
									   index=columns, columns=index).transpose()

	return lists_frame


def writeDFsToExcelSheets(excel_name, data_frames, sheet_names,
                          read_files=False, index=True):
	""" Takes in a list of dataframe fileNames, \
	and writes these to an excel sheet. Automatically adds an extra column to \
	DESingle output called Type_renamed, making the following conversion from \
	the Type column: typeToType2 = {'DEs': 'DEa', 'DEa': 'DEm', 'DEg': 'DE',
														   numpy.nan: numpy.nan}

	Args:
		excel_name (str): Specifies the name of the excel sheet.

		data_frames (list-like<str> or list-like(pandas.DataFrame)):
									List of strings specifying locations of \
									dataframes to read if read_files=True; \
									otherwise list of dataframes to write \
									into excel sheets.

		sheet_names (list-like<str>): List of strings specifying the names of \
															 excel sheet to add.

		read_files (bool): Whether to read in the file names as pandas data \
						frames, or these are already pandas data frames.
	"""

	# Create a Pandas Excel writer using XlsxWriter as the engine.
	writer = pandas.ExcelWriter(excel_name,
								engine='xlsxwriter')
	workbook = writer.book

	for i, data_frame in enumerate(data_frames):
		if read_files:
			df = pandas.read_csv(data_frame, sep='\t', header=0)

		else:
			df = data_frame

		sheetName = sheet_names[i]

		# Convert the dataframe to an XlsxWriter Excel object.
		df.to_excel(writer, sheet_name=sheetName, startrow=0,
					index=index)

		# Extracting the sheet for extra formatting
		sheet = writer.sheets[sheetName]

	# Close the Pandas Excel writer and output the Excel file.
	writer.save()


