# Molecule Align
Align molecules using the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm).

## Description

This tool aligns atoms across different snapshots in a trajectory using the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm), allowing only rotational adjustments without any stretching or shearing.

Users can select either all atoms or a subset as the reference for alignment. Optionally, a reference frame can be specified, defaulting to the first frame of the trajectory. This frame serves as the baseline orientation for the atoms, with their [centroid](https://en.wikipedia.org/wiki/Centroid) repositioned to the origin.

### Requirements

This modifier requires unique identifiers for each atom. If these are not available, they will be generated during runtime. In that case the modifier assumes a consistent atom order throughout the trajectory.

It's crucial to provide an unwrapped trajectory to ensure accurate behavior. If your trajectory is wrapped, use the [unwrap trajectories modifier](https://www.ovito.org/docs/current/reference/pipelines/modifiers/unwrap_trajectories.html) before applying this modifier. The behavior with wrapped trajectories is unpredictable.

When operating on a subset of atoms, ensure the selection remains consistent throughout the trajectory. Utilize the [freeze property](https://www.ovito.org/docs/current/reference/pipelines/modifiers/freeze_property.html#particles-modifiers-freeze-property) modifier to maintain a consistent selection, or the modifier's output may be undefined.

### Output

The modifier calculates the root mean square deviation (RMSD) of atomic coordinates post-alignment relative to the reference frame. The RMSD for selected atoms is recorded in the `MoleculeAlign.RMSD` [global attribute](https://www.ovito.org/docs/current/reference/data_inspector/attributes.html). The RMSD for all atoms is available in the `MoleculeAlign.RMSD_all` attribute.

Additionally, it provides an RMSD value for each atom and generates a data table showing RMSD values across frames, which populates as the trajectory is analyzed in OVITO Pro.

## Parameters 

- `only_selected` / "Use only selected particles": This option allows alignment using only selected particles as the reference. `Default = True`
- `reference_frame` / "Reference frame": Specifies the frame to use as the reference for alignment. `Default = 0`

## Example

Here's an example illustrating how the modifier aligns the highlighted (green) atoms throughout the trajectory.

![Example image showing the molecule alignment](examples/Example_01.png)

![Example video showing the molecule alignment](examples/Example_01.mp4)

## Installation
- OVITO Pro [integrated Python interpreter](https://docs.ovito.org/python/introduction/installation.html#ovito-pro-integrated-interpreter):
  ```
  ovitos -m pip install --user git+https://github.com/ovito-org/MoleculeAlign.git
  ``` 
  The `--user` option is recommended and [installs the package in the user's site directory](https://pip.pypa.io/en/stable/user_guide/#user-installs).

- Other Python interpreters or Conda environments:
  ```
  pip install git+https://github.com/ovito-org/MoleculeAlign.git
  ```

## Technical information / dependencies
- Tested on OVITO version 3.10.3

## Contact
Daniel Utt (utt@ovito.org)
