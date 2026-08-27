# WALLABY source-only TAP access links

Post-Stage-0 acquisition helper. These links query **only** the public 30-arcsec WALLABY DR2 source catalogue. No kinematic catalogue is referenced.

- [Source catalogue schema](https://ws-uv.canfar.net/youcat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+column_name%2Cdatatype%2Cdescription%2Cunit+FROM+TAP_SCHEMA.columns+WHERE+table_name%3D%27cirada.Wallaby_dr2_source_catalogue%27)
- [Five-row source-catalogue sample](https://ws-uv.canfar.net/youcat/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=SELECT+TOP+5+%2A+FROM+cirada.Wallaby_dr2_source_catalogue)

These links were added only after the Stage-0 preregistration commit `15abeac9b6d285ce42caf26498115c0686d9bcd3` and do not alter that earlier timestamp.
