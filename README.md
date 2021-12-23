# WAVEPACK

`WAVEPACK`: routines to compute the wavelet transform of a time series, and significance levels.

## Get Started

### Dependencies

- Git
- [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

Fortran Package Manager (fpm) is a package manager and build system for Fortran.  
You can build `wavepack` using provided `fpm.toml`:

```sh
fpm build --profile release
```

To use `wavepack` within your `fpm` project, add the following to your `fpm.toml` file:

```toml
[dependencies]
wavepack = { git="https://github.com/zoziha/wavepack" }
```

## Links

- [danielrios12/wavelet---torrence-e-compo](https://github.com/danielrios12/wavelet---torrence-e-compo)
- [fortran-lang/fftpack](https://github.com/fortran-lang/fftpack)
- [netlib/random](http://www.netlib.org/random/)