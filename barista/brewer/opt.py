try: 
    from barista.brewer.BrewerJobManager import BrewerJobManager
except:
    from BrewerJobManager import BrewerJobManager 

from contextlib import redirect_stdout

job = BrewerJobManager('brewer.TEMPLATE')

with open('brewer_opt.log', 'w') as logfile:
    with redirect_stdout(logfile):
        job.run_job()
