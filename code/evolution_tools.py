#----------------------------------------------------------------------------------------
# Modules for processing evolutionary data.
#----------------------------------------------------------------------------------------

import io
import re
import numpy as np

def read_consurf_score_file (inPath, col = None):
    """Read ConSurf score file.

    Args:
        inPath (Path): path to file containing ConSurf scores.
        col (list): list of columns to read, other than POS column.
                    Defaults are columns SEQ and SCORE.

    """
    if not col:
        col = ['SEQ', 'SCORE']
    allcol = ['POS', 'SEQ', 'SCORE', 'COLOR', 'CONFIDENCE INTERVAL', 'CONFIDENCE INTERVAL COLORS',
                'B/E', 'FUNCTION', 'MSA DATA', 'RESIDUE VARIETY']
    
    scores = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        for i in range(16):
            next(f)
        for line in f:
            line = line.strip()
            if line:
                ls = re.split('\t+', line)
                ls = [c.replace(' ', '') for c in ls]
                val = {c:v for c, v in zip(allcol[1:], ls[1:])}
                try:
                    val['SCORE'] = float(val['SCORE'])
                except:
                    val['SCORE'] = np.nan
                scores[int(ls[0])] = {k:v for k,v in val.items() if k in col}
            else:
                break
    return scores
