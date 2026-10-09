#!/usr/bin/env bash
# The model environment of the published pipeline (AMR_benchmarking install/amr_env, as built for
# the Hu et al. 2024 reproduction): conda python 3.7 + pip-installed scikit-learn 0.24.1,
# numpy 1.20.1, pandas 1.2.2, scipy 1.6.0. Plus PyYAML and s2g2p (no deps), which is all
# s2g2p needs to run `evaluate --preset hu2024` under the published library versions.
#   bash containers/hu2024ml/make-env.sh <prefix>
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
prefix=${1:?prefix}
micromamba create -y -p "$prefix" --file "$here/amr_env.explicit.txt" 2>/dev/null \
  || conda create -y -p "$prefix" --file "$here/amr_env.explicit.txt"
"$prefix/bin/pip" install --no-deps -r "$here/amr_env.pip.txt"
"$prefix/bin/pip" install --no-deps "pyyaml==6.0.1"
# s2g2p by path, not by pip: its build backend (hatchling) needs Python >= 3.8.
site=$("$prefix/bin/python" -c "import site; print(site.getsitepackages()[0])")
echo "$(cd "$here/../.." && pwd)/src" > "$site/s2g2p.pth"
"$prefix/bin/python" -c "import sklearn, numpy, pandas; assert sklearn.__version__ == '0.24.1', sklearn.__version__; print('ok', sklearn.__version__, numpy.__version__, pandas.__version__)"
