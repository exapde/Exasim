class Exasim < Formula
  desc "Exasim PDE simulation framework"
  homepage "https://github.com/exapde/Exasim"
  url "https://github.com/exapde/Exasim/archive/refs/tags/v1.4.tar.gz"
  sha256 "REPLACE_WITH_V1_4_TARBALL_SHA256"
  license "MIT"

  depends_on "cmake" => :build
  depends_on "open-mpi"
  depends_on "openblas"

  def install
    ENV.prepend_path "PATH", Formula["open-mpi"].opt_bin

    # Build Kokkos libraries
    Dir.chdir("kokkos") do
      system "make", "-f", "Makefile.builds", "serial"
    end

    # Build Metis and ParMetis libraries
    Dir.chdir("metis") do
      system "make", "metis"
    end

    # Build Text2Code executable
    Dir.chdir("text2code") do
      system "make", "text2code"
    end

    # Build Exasim executables
    mkdir_p "build"
    Dir.chdir("build") do
      system "cmake",
             "-D", "EXASIM_NOMPI=ON",
             "-D", "EXASIM_MPI=ON",
             "-D", "EXASIM_CUDA=OFF",
             "-D", "EXASIM_HIP=OFF",
             "-D", "WITH_TEXT2CODE=ON",
             "-D", "WITH_PARMETIS=ON",
             "../install"
      system "cmake", "--build", "."
    end

    # Install generated executables if present.
    bin.install "build/text2code" if (buildpath/"build/text2code").exist?
    bin.install "build/t2cusecmake" if (buildpath/"build/t2cusecmake").exist?
    bin.install "build/cput2cEXASIM" if (buildpath/"build/cput2cEXASIM").exist?
    bin.install "build/cpumpit2cEXASIM" if (buildpath/"build/cpumpit2cEXASIM").exist?
  end

  test do
    exe = bin/"cput2cEXASIM"
    assert_predicate exe, :exist?
    out = shell_output("#{exe} 2>&1", 1)
    assert_match "Usage: ./Exasim <pdeapp.txt>", out
  end
end
