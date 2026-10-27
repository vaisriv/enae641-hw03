{ pkgs, inputs, ... }:
inputs.treefmt-nix.lib.mkWrapper pkgs {
    projectRootFile = "flake.nix";

    # nix
    programs = {
        deadnix.enable = true;
        nixfmt = {
            enable = true;
            indent = 4;
            strict = true;
        };
    };

    # md/yaml
    programs.prettier = {
        enable = true;
        settings = {
            tabWidth = 4;
        };
    };

    # typst
    programs.typstyle = {
        enable = true;
        indentWidth = 4;
    };
}
