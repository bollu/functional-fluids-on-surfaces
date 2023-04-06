{ pkgs ? import <nixpkgs> {} }:
  let
    my-python-packages = ps: with ps; [
      pandas
      requests
      # other python packages
    ];
  in

  pkgs.mkShell {
    nativeBuildInputs = with pkgs; [
      gcc
      xlibsWrapper
      xorg.libX11
      xorg.libXrandr
      xorg.libXcursor
      xorg.libXi
      stdenv.cc.cc.lib
      zlib
    ];
    packages = [
        (pkgs.python3.withPackages my-python-packages) # we have defined this in the installation section
    ];
    buildInputs = with pkgs; [
      gcc
      clang
      stdenv.cc.cc.lib
      xlibsWrapper
      xorg.libX11
      xorg.libXrandr
      xorg.libXcursor
      xorg.libXi
      glibc.static
      zlib
    ];
    shellHook = ''
      export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${pkgs.lib.makeLibraryPath [pkgs.xorg.libX11 pkgs.xlibsWrapper pkgs.stdenv.cc.cc.lib pkgs.zlib]}"
    '';

}

