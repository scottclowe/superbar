CHANGELOG
=========

## v1.5.0 - 2016-10-16

### Added
- PStarLatex option allows changing whether significance text is rendered
  using the LaTeX interpreter, with the default now as `'off'`.
- PStarIcon option allows changing the symbol used for stars.


## v1.4.1 - 2016-10-16

### Added
- Support for R2007a-R2013a `inputParser` method.

### Fixed
- Avoid using missing `narginchk` in pre-2011b.


## v1.4.0 - 2016-07-24

### Added
- Support for cell array colour inputs.
- Support inheriting errorbar colour from char string with more than one
  element.
- BarLineWidth input added, with default of 2.

### Changed
- Bar faces now default to transparent if a bar edge colour is given.
- Default errorbar colour inherited from bar face colour will now be lighter
  than the bar if the sum of the face RGB values is less than 0.65.


## v1.3.2 - 2016-07-24
### Added
- Tutorial for superbar API.
- README.md.
- Logo.

### Changed
- Default PLineSourceRelativeSpacing now based on errorbar width, if present.
- Default PStarThreshold is now `[0.05, 0.01, 0.001, 0.0001]`.

### Fixed
- Behaviour of > sign when showing significance. Now > appears instead of a
  star when all thresholds are exceeded.


## v1.3.1 - 2016-07-24

### Added
- Add unit tests for SUPERBAR.

### Fixed
- Fix bug when plotting into non-current axes.


## v1.3.0 - 2016-07-23

### Added
- Add unit tests for supererr.

### Changed
- When there is no source separation, offset line starts.


## v1.2.1 - 2016-06-27

### Changed
- Change output shapes of HB, HPT, HPL, HPB to match the shape of Y.


## v1.2.0 - 2016-06-14

### Added
- Show "n.s." when not significant
- Background colour to P-value text

### Changed
- P-value text in LaTeX mode
- Change shape of outputs hpl, hpt, hpb
- Update defaults for p thresholds, errorbar widths, fontsize

### Fixed
- Avoid axes input to text and line fns, supporting pre-R2016a (tested in R2015b)


## v1.1.1 - 2016-06-14

### Changed
- Change default orientation of star text for horizontal bars; now normal for
single p-values and sideways for pair-wise p-values.


## v1.1.0 - 2016-06-10

### Added
- Show pairwise comparison lines on horizontally oriented plots.


## v1.0.1 - 2016-06-05

### Fixed
- Bug fixes.


## v1.0.0 - 2016-06-03
- Initial release.
