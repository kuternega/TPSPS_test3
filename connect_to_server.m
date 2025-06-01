
f = parfeval(@() connect(), 0);
while strcmp(f.State, 'running') || strcmp(f.State, 'queued')
    pause(0.01);
end
connect()

function connect()
   addpath('Z:\COMSOL56\Multiphysics\mli'); % путь к mli дирректории COMSOL
   mphstart('localhost', 2036);
end

