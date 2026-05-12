function results = exasimtest(np)
%EXASIMTEST Cross-frontend and backend ABI smoke tests for Exasim.
%
% RESULTS = EXASIMTEST() exercises one Poisson problem through:
%   1. Matlab, Julia, and Python frontends in examples/Poisson/poisson2d.
%   2. Text2Code and BuiltIn dynamic-library ABI executables for
%      apps/poisson/poisson2d.
%   3. Header-only Kokkos-kernel ABI executables for apps/poisson/poisson2d.
%
% Optional environment overrides:
%   EXASIM_CMAKE   path to cmake
%   EXASIM_MPIRUN  path to mpirun
%   EXASIM_TEXT2CODE path to text2code
%   EXASIM_JULIA   path to julia
%   EXASIM_PYTHON  path to python3
%   EXASIM_NP      MPI process count for backend runs

testFile = mfilename('fullpath');
testDir = fileparts(testFile);
exasimDir = fileparts(testDir);

exampleDir = fullfile(exasimDir, 'examples', 'Poisson', 'poisson2d');
appDir = fullfile(exasimDir, 'apps', 'poisson', 'poisson2d');
buildDir = fullfile(exasimDir, 'build');
installDir = fullfile(exasimDir, 'install');
kokkosKernelHeader = fullfile(appDir, 'poisson2d.hpp');

cmake = getenvDefault('EXASIM_CMAKE', findExecutable('cmake', '/opt/homebrew/bin/cmake'));
mpirun = getenvDefault('EXASIM_MPIRUN', findExecutable('mpirun', '/opt/homebrew/bin/mpirun'));
text2code = getenvDefault('EXASIM_TEXT2CODE', findExecutable('text2code', fullfile(buildDir, 'text2code')));
julia = getenvDefault('EXASIM_JULIA', findJuliaExecutable());
python = getenvDefault('EXASIM_PYTHON', 'python3');
if isnan(np) || np < 1
    error('np must be a positive integer.');
end
np = round(np);

results = struct('name', {}, 'command', {}, 'status', {});

fprintf('==> Exasim smoke tests\n');
fprintf('Exasim root: %s\n', exasimDir);
fprintf('Build path : %s\n', buildDir);
fprintf('MPI ranks  : %d\n\n', np);

%% Frontend tests
resetCMakeCache(buildDir);
results(end+1) = runMatlabScript('frontend-matlab-poisson2d', exampleDir, 'pdeapp.m'); %#ok<SAGROW>
resetCMakeCache(buildDir);

results(end+1) = runShellTest('frontend-julia-poisson2d', exampleDir, ...
    sprintf('%s %s', shellQuote(julia), shellQuote(fullfile(exampleDir, 'pdeapp.jl')))); %#ok<SAGROW>

resetCMakeCache(buildDir);
results(end+1) = runShellTest('frontend-python-poisson2d', exampleDir, ...
    sprintf('%s %s', shellQuote(python), shellQuote(fullfile(exampleDir, 'pdeapp.py')))); %#ok<SAGROW>

%% Backend ABI build
if ~exist(buildDir, 'dir')
    mkdir(buildDir);
end
resetCMakeCache(buildDir);

text2codeGenerate = sprintf('%s %s', shellQuote(text2code), shellQuote('pdeapp.txt'));
results(end+1) = runShellTest('backend-text2code-generate', appDir, text2codeGenerate); %#ok<SAGROW>

cmakeConfigure = strjoin({ ...
    shellQuote(cmake), ...
    '-S', shellQuote(installDir), ...
    '-B', shellQuote(buildDir), ...
    '-D', 'EXASIM_NOMPI=ON', ...
    '-D', 'EXASIM_MPI=ON', ...
    '-D', 'WITH_PARMETIS=ON', ...
    '-D', 'WITH_TEXT2CODE=ON', ...
    '-D', 'WITH_BUILTINMODEL=ON', ...
    '-D', ['KOKKOSKERNEL_HEADER=' shellQuote(kokkosKernelHeader)]}, ' ');

results(end+1) = runShellTest('backend-cmake-configure', exasimDir, cmakeConfigure); %#ok<SAGROW>

backendTargets = { ...
    'cput2cEXASIM', ...
    'cpumpit2cEXASIM', ...
    'cpubuiltinEXASIM', ...
    'cpumpibuiltinEXASIM', ...
    'cpukkEXASIM', ...
    'cpumpikkEXASIM'};

cmakeBuild = sprintf('%s --build %s --target %s', ...
    shellQuote(cmake), shellQuote(buildDir), strjoin(backendTargets, ' '));
results(end+1) = runShellTest('backend-cmake-build', exasimDir, cmakeBuild); %#ok<SAGROW>

%% Backend ABI run tests
results(end+1) = runBackendExecutable('backend-cput2c', appDir, buildDir, 'cput2cEXASIM', '', 'pdeapp.txt'); %#ok<SAGROW>
results(end+1) = runBackendExecutable('backend-cpumpit2c', appDir, buildDir, 'cpumpit2cEXASIM', mpirun, 'pdeapp.txt', np); %#ok<SAGROW>
results(end+1) = runBackendExecutable('backend-cpubuiltin', appDir, buildDir, 'cpubuiltinEXASIM', '', 'pdeapp.txt'); %#ok<SAGROW>
results(end+1) = runBackendExecutable('backend-cpumpibuiltin', appDir, buildDir, 'cpumpibuiltinEXASIM', mpirun, 'pdeapp.txt', np); %#ok<SAGROW>
results(end+1) = runBackendExecutable('backend-cpukk', appDir, buildDir, 'cpukkEXASIM', '', 'pdeapp.txt'); %#ok<SAGROW>
results(end+1) = runBackendExecutable('backend-cpumpikk', appDir, buildDir, 'cpumpikkEXASIM', mpirun, 'pdeapp.txt', np); %#ok<SAGROW>

%% Summary
fprintf('\n==> Exasim test summary\n');
failed = false;
for i = 1:numel(results)
    if results(i).status == 0
        mark = 'PASS';
    else
        mark = 'FAIL';
        failed = true;
    end
    fprintf('%-6s %s\n', mark, results(i).name);
end

if failed
    error('One or more Exasim smoke tests failed.');
end

fprintf('All Exasim smoke tests passed.\n');
end

function result = runMatlabScript(name, workDir, scriptName)
    fprintf('\n==> %s\n', name);
    fprintf('cd %s\n', workDir);
    fprintf('run %s\n', scriptName);

    oldDir = pwd();
    cleanup = onCleanup(@() cd(oldDir));
    result = struct('name', name, 'command', ['run ' scriptName], 'status', 1);
    try
        cd(workDir);
        run(scriptName);
        result.status = 0;
    catch err
        fprintf(2, 'FAILED: %s\n', name);
        fprintf(2, '%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
        result.status = 1;
    end
end

function result = runShellTest(name, workDir, command)
    fprintf('\n==> %s\n', name);
    fprintf('cd %s\n', workDir);
    fprintf('%s\n', command);

    oldDir = pwd();
    cleanup = onCleanup(@() cd(oldDir));
    cd(workDir);
    status = system(command);
    result = struct('name', name, 'command', command, 'status', status);
end

function result = runBackendExecutable(name, workDir, buildDir, exeName, mpirun, inputFile, np)
    if nargin < 7
        np = 1;
    end

    exePath = fullfile(buildDir, exeName);
    if isempty(mpirun)
        command = sprintf('%s %s', shellQuote(exePath), shellQuote(inputFile));
    else
        command = sprintf('%s -np %d %s %s', shellQuote(mpirun), np, shellQuote(exePath), shellQuote(inputFile));
    end

    result = runShellTest(name, workDir, command);
end

function value = getenvDefault(name, defaultValue)
    value = getenv(name);
    if isempty(value)
        value = defaultValue;
    end
end

function exe = findExecutable(name, fallback)
    [status, output] = system(sprintf('command -v %s', shellQuote(name)));
    if status == 0
        exe = strtrim(output);
    elseif exist(fallback, 'file')
        exe = fallback;
    else
        exe = name;
    end
end

function exe = findJuliaExecutable()
    [status, output] = system('command -v julia');
    if status == 0
        exe = strtrim(output);
        return;
    end

    candidates = {};
    homeDir = getenv('HOME');
    if ~isempty(homeDir)
        candidates{end+1} = fullfile(homeDir, '.juliaup', 'bin', 'julia'); %#ok<AGROW>
        candidates{end+1} = fullfile(homeDir, 'bin', 'julia'); %#ok<AGROW>
    end
    candidates{end+1} = '/opt/homebrew/bin/julia';
    candidates{end+1} = '/usr/local/bin/julia';

    appCandidates = dir('/Applications/Julia*.app/Contents/Resources/julia/bin/julia');
    for i = 1:numel(appCandidates)
        candidates{end+1} = fullfile(appCandidates(i).folder, appCandidates(i).name); %#ok<AGROW>
    end

    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file')
            exe = candidates{i};
            return;
        end
    end

    exe = 'julia';
end

function q = shellQuote(value)
    value = char(value);
    q = ['''' strrep(value, '''', '''"''"''') ''''];
end

function resetCMakeCache(buildDir)
    cacheFile = fullfile(buildDir, 'CMakeCache.txt');
    cmakeFilesDir = fullfile(buildDir, 'CMakeFiles');

    if exist(cacheFile, 'file')
        delete(cacheFile);
    end
    if exist(cmakeFilesDir, 'dir')
        rmdir(cmakeFilesDir, 's');
    end
end
