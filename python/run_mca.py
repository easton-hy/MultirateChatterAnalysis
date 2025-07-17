import datetime
import os
from .multirate_analysis_cs import multirate_analysis
from .save_all import save_results


def run():
    supertime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    dirname = os.path.join('result', f'SLD_MCA_{supertime}')
    tool_list = [1, 2]
    for tool in tool_list:
        sdm = multirate_analysis(tool)
        save_results(dirname, tool, sdm)

if __name__ == '__main__':
    run()
