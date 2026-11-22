# Split statistics (Role 1)

Generated 2026-11-20 by `scripts/make_splits.py`. Residue counts use `valid_mask` (ordered residues with a Q8 label; on CB513 also the evaluation mask).

- CB513 leakage rule (Train proteins removed): E-value <= 0.001, identity >= 25%, alignment >= 20 residues.
- Train grouping rule (same fold): E-value <= 0.001, identity >= 30%, alignment >= 20 residues, covering >= 50% of the shorter protein.
- Status: **complete (exact/substring + MMseqs2 CB513-vs-Train and Train-vs-Train)**; 9,376 groups, largest 713 proteins.

| Partition | Proteins | Valid residues | Median len | Max len | >1022 | G | G % | H | H % | G next to H (% of G) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fold 0 | 2,124 | 450,321 | 208 | 1009 | 0 | 17,167 | 3.81 | 156,163 | 34.68 | 3,605 (21.0%) |
| fold 1 | 1,831 | 450,321 | 224 | 1632 | 8 | 17,171 | 3.81 | 155,211 | 34.47 | 3,656 (21.3%) |
| fold 2 | 1,842 | 450,321 | 216 | 1475 | 8 | 17,630 | 3.91 | 153,085 | 33.99 | 3,640 (20.6%) |
| fold 3 | 1,836 | 450,321 | 221 | 1267 | 7 | 17,401 | 3.86 | 152,307 | 33.82 | 3,679 (21.1%) |
| fold 4 | 1,838 | 450,320 | 222 | 1491 | 11 | 17,036 | 3.78 | 156,906 | 34.84 | 3,511 (20.6%) |
| all CV folds | 9,471 | 2,251,604 | 217 | 1632 | 34 | 86,405 | 3.84 | 773,672 | 34.36 | 18,091 (20.9%) |
| validation (-1) | 1,053 | 245,405 | 211 | 1108 | 1 | 9,484 | 3.86 | 85,588 | 34.88 | 2,056 (21.7%) |
| CB513 test (-2) | 513 | 83,772 | 238 | 874 | 0 | 3,098 | 3.70 | 25,706 | 30.69 | 555 (17.9%) |
| excluded | 324 | 101,838 | 294 | 1073 | 1 | 4,650 | 4.57 | 32,302 | 31.72 | 761 (16.4%) |

## Q8 distribution (% of valid residues)

| Partition | G | H | I | B | E | S | T | C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| fold 0 | 3.81 | 34.68 | 0.54 | 0.93 | 21.90 | 7.85 | 10.92 | 19.36 |
| fold 1 | 3.81 | 34.47 | 0.58 | 1.05 | 21.59 | 8.15 | 10.84 | 19.51 |
| fold 2 | 3.91 | 33.99 | 0.55 | 1.03 | 21.77 | 8.23 | 10.79 | 19.71 |
| fold 3 | 3.86 | 33.82 | 0.51 | 1.08 | 21.82 | 8.16 | 11.09 | 19.64 |
| fold 4 | 3.78 | 34.84 | 0.55 | 1.01 | 21.40 | 8.17 | 10.86 | 19.39 |
| all CV folds | 3.84 | 34.36 | 0.55 | 1.02 | 21.70 | 8.11 | 10.90 | 19.52 |
| validation (-1) | 3.86 | 34.88 | 0.58 | 0.94 | 21.76 | 7.95 | 10.82 | 19.21 |
| CB513 test (-2) | 3.70 | 30.69 | 0.57 | 1.32 | 21.05 | 9.83 | 11.65 | 21.21 |
| excluded | 4.57 | 31.72 | 0.49 | 1.23 | 21.51 | 8.70 | 12.01 | 19.77 |

## CB513 vs Train similarity (MMseqs2, E-value <= 1e-3, alignment >= 20)

NetSurfP-2.0 states that no training protein has > 25% identity to CB513. MMseqs2 local alignments still find the matches below, so these Train proteins are excluded.

| Identity >= | Pairs | CB513 chains | Train chains |
|---:|---:|---:|---:|
| 25% | 470 | 203 | 323 |
| 30% | 202 | 108 | 154 |
| 40% | 96 | 47 | 66 |
| 50% | 72 | 27 | 45 |
| 90% | 42 | 13 | 29 |

## Excluded training proteins

| uid | pdb_id | reason |
|---|---|---|
| train_00046 | 1bf2-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.25, alignment 213, E=4.2e-06 |
| train_00054 | 1brt-A | homolog of CB513 cb513_00106 (1edd-A): identity 0.34, alignment 85, E=1.5e-04 |
| train_00065 | 1c3c-A | homolog of CB513 cb513_00129 (1fuq-B): identity 0.30, alignment 198, E=1.2e-11 |
| train_00084 | 1cs0-G | homolog of CB513 cb513_00031 (1bnc-B): identity 0.26, alignment 174, E=4.8e-13 |
| train_00090 | 1cyw-A | homolog of CB513 cb513_00076 (1cyx-A): identity 0.97, alignment 205, E=4.8e-132 |
| train_00103 | 1d4t-A | homolog of CB513 cb513_00018 (1aya-B): identity 0.31, alignment 96, E=7.3e-07 |
| train_00147 | 1ed1-A | homolog of CB513 cb513_00166 (1hiw-S): identity 0.60, alignment 100, E=2.1e-34 |
| train_00158 | 1ehy-B | homolog of CB513 cb513_00106 (1edd-A): identity 0.33, alignment 96, E=3.8e-06 |
| train_00244 | 1g5a-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.29, alignment 180, E=9.6e-09 |
| train_00261 | 1gk6-A | homolog of CB513 cb513_00136 (1gcm-C): identity 0.72, alignment 32, E=1.1e-11 |
| train_00305 | 1h2b-A | homolog of CB513 cb513_00507 (8adh-A): identity 0.26, alignment 367, E=3.7e-21 |
| train_00439 | 1jfb-A | homolog of CB513 cb513_00496 (6cpp-A): identity 0.30, alignment 293, E=1.1e-27 |
| train_00443 | 1jhc-A | homolog of CB513 cb513_00347 (1umu-B): identity 0.32, alignment 113, E=2.2e-12 |
| train_00453 | 1jm1-A | homolog of CB513 cb513_00293 (1rie-A): identity 0.36, alignment 74, E=4.5e-06 |
| train_00499 | 1k8k-C | homolog of CB513 cb513_00152 (1gp2-B): identity 0.26, alignment 121, E=7.6e-05 |
| train_00504 | 1k8u-A | homolog of CB513 cb513_00462 (3icb-A): identity 0.45, alignment 62, E=1.1e-08 |
| train_00506 | 1ka1-A | homolog of CB513 cb513_00407 (2hhm-B): identity 0.27, alignment 194, E=9.6e-07 |
| train_00525 | 1kol-A | homolog of CB513 cb513_00507 (8adh-A): identity 0.28, alignment 266, E=1.3e-17 |
| train_00563 | 1lf7-A | homolog of CB513 cb513_00285 (1rbp-A): identity 0.27, alignment 113, E=3.7e-04 |
| train_00614 | 1mgt-A | homolog of CB513 cb513_00309 (1sfe-A): identity 0.37, alignment 77, E=5.6e-06 |
| train_00626 | 1ms9-A | homolog of CB513 cb513_00116 (1euu-A): identity 0.25, alignment 231, E=5.6e-08 |
| train_00631 | 1mtz-A | homolog of CB513 cb513_00106 (1edd-A): identity 0.27, alignment 107, E=8.2e-04 |
| train_00632 | 1mun-A | homolog of CB513 cb513_00371 (2abk-A): identity 0.30, alignment 127, E=2.7e-10 |
| train_00725 | 1o7f-A | homolog of CB513 cb513_00132 (1g6n-A): identity 0.28, alignment 81, E=9.1e-05 |
| train_00827 | 1psr-B | homolog of CB513 cb513_00462 (3icb-A): identity 0.31, alignment 68, E=1.1e-04 |
| train_00830 | 1puc-A | homolog of CB513 cb513_00055 (1cks-C): identity 0.52, alignment 81, E=5.7e-23 |
| train_00889 | 1qq5-A | homolog of CB513 cb513_00193 (1jud-A): identity 0.48, alignment 193, E=1.9e-52 |
| train_00902 | 1qwd-B | homolog of CB513 cb513_00021 (1bbp-A): identity 0.28, alignment 132, E=8.1e-07 |
| train_00913 | 1r0m-D | homolog of CB513 cb513_00226 (1mns-A): identity 0.27, alignment 260, E=6.3e-15 |
| train_00936 | 1req-B | homolog of CB513 cb513_00289 (1req-C): identity 0.30, alignment 314, E=1.2e-31 |
| train_00941 | 1ri5-A | homolog of CB513 cb513_00362 (1xva-B): identity 0.28, alignment 156, E=7.4e-04 |
| train_00978 | 1ryo-A | homolog of CB513 cb513_00249 (1ovb-A): identity 0.56, alignment 156, E=1.4e-48 |
| train_01079 | 1tii-C | identical sequence to cb513_00331 (1tii-C) |
| train_01083 | 1tkj-A | homolog of CB513 cb513_00008 (1amp-A): identity 0.29, alignment 173, E=1.9e-10 |
| train_01091 | 1tu5-A | homolog of CB513 cb513_00245 (1oac-B): identity 0.27, alignment 427, E=7.3e-33 |
| train_01155 | 1uh4-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.25, alignment 480, E=2.0e-17 |
| train_01228 | 1vg0-A | homolog of CB513 cb513_00146 (1gnd-A): identity 0.28, alignment 287, E=4.0e-28 |
| train_01259 | 1vl1-A | homolog of CB513 cb513_00172 (1hor-B): identity 0.29, alignment 114, E=1.2e-06 |
| train_01304 | 1w07-A | homolog of CB513 cb513_00465 (3mdd-B): identity 0.29, alignment 220, E=9.9e-12 |
| train_01321 | 1w7c-A | homolog of CB513 cb513_00245 (1oac-B): identity 0.26, alignment 501, E=1.3e-33 |
| train_01324 | 1w98-B | homolog of CB513 cb513_00123 (1fin-D): identity 0.32, alignment 196, E=1.8e-26 |
| train_01338 | 1wfb-A | identical sequence to cb513_00360 (1wfb-B) |
| train_01343 | 1wkr-A | homolog of CB513 cb513_00490 (5er2-E): identity 0.31, alignment 233, E=4.8e-15 |
| train_01344 | 1wle-B | homolog of CB513 cb513_00307 (1ses-A): identity 0.34, alignment 432, E=7.7e-59 |
| train_01386 | 1wzn-A | homolog of CB513 cb513_00362 (1xva-B): identity 0.27, alignment 115, E=2.0e-08 |
| train_01400 | 1x9m-A | homolog of CB513 cb513_00200 (1ktq-A): identity 0.26, alignment 367, E=7.8e-19 |
| train_01459 | 1xu9-D | homolog of CB513 cb513_00121 (1fdt-A): identity 0.28, alignment 165, E=1.4e-06 |
| train_01511 | 1yht-A | homolog of CB513 cb513_00281 (1qbb-A): identity 0.27, alignment 199, E=9.0e-11 |
| train_01559 | 1z0x-B | homolog of CB513 cb513_00444 (2trt-A): identity 0.27, alignment 148, E=6.7e-11 |
| train_01578 | 1z9l-A | homolog of CB513 cb513_00230 (1msp-B): identity 0.28, alignment 80, E=1.8e-07 |
| train_01726 | 2b9c-A | homolog of CB513 cb513_00136 (1gcm-C): identity 0.68, alignment 31, E=2.8e-10 |
| train_01738 | 2bf6-A | homolog of CB513 cb513_00116 (1euu-A): identity 0.28, alignment 436, E=3.0e-28 |
| train_01755 | 2bl9-A | homolog of CB513 cb513_00498 (6dfr-A): identity 0.30, alignment 104, E=9.5e-10 |
| train_01775 | 2btj-D | homolog of CB513 cb513_00139 (1gfl-B): identity 0.31, alignment 139, E=4.7e-20 |
| train_01811 | 2cdc-B | homolog of CB513 cb513_00507 (8adh-A): identity 0.32, alignment 97, E=4.5e-08 |
| train_02056 | 2fba-A | homolog of CB513 cb513_00018 (1aya-B): identity 0.34, alignment 59, E=7.4e-04 |
| train_02069 | 2fdb-N | homolog of CB513 cb513_00026 (1bfg-A): identity 0.30, alignment 121, E=4.3e-17 |
| train_02074 | 2fel-B | homolog of CB513 cb513_00129 (1fuq-B): identity 0.27, alignment 205, E=1.2e-04 |
| train_02094 | 2fj8-A | homolog of CB513 cb513_00321 (1tab-I): identity 0.43, alignment 46, E=8.4e-06 |
| train_02116 | 2fre-A | homolog of CB513 cb513_00241 (1nox-A): identity 0.30, alignment 187, E=1.5e-04 |
| train_02175 | 2gdq-B | homolog of CB513 cb513_00226 (1mns-A): identity 0.28, alignment 213, E=1.0e-18 |
| train_02215 | 2grv-B | homolog of CB513 cb513_00420 (2olb-A): identity 0.25, alignment 151, E=1.3e-06 |
| train_02286 | 2him-A | homolog of CB513 cb513_00458 (3eca-B): identity 0.28, alignment 269, E=2.6e-14 |
| train_02328 | 2hxi-B | homolog of CB513 cb513_00444 (2trt-A): identity 0.36, alignment 85, E=9.2e-07 |
| train_02329 | 2hxo-A | homolog of CB513 cb513_00444 (2trt-A): identity 0.27, alignment 155, E=6.0e-09 |
| train_02331 | 2hxt-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.31, alignment 220, E=2.6e-17 |
| train_02395 | 2id4-A | homolog of CB513 cb513_00301 (1s01-A): identity 0.30, alignment 265, E=5.3e-09 |
| train_02434 | 2irp-B | homolog of CB513 cb513_00127 (1fua-A): identity 0.33, alignment 99, E=7.1e-10 |
| train_02436 | 2isk-G | homolog of CB513 cb513_00241 (1nox-A): identity 0.25, alignment 190, E=6.4e-05 |
| train_02480 | 2j91-B | homolog of CB513 cb513_00129 (1fuq-B): identity 0.26, alignment 181, E=2.9e-04 |
| train_02482 | 2j9o-C | homolog of CB513 cb513_00451 (3blm-A): identity 0.26, alignment 247, E=1.1e-07 |
| train_02513 | 2jhj-A | homolog of CB513 cb513_00371 (2abk-A): identity 0.30, alignment 126, E=5.3e-04 |
| train_02514 | 2jhn-A | homolog of CB513 cb513_00371 (2abk-A): identity 0.30, alignment 126, E=5.3e-04 |
| train_02520 | 2jku-A | homolog of CB513 cb513_00023 (1bdo-A): identity 0.36, alignment 72, E=9.7e-09 |
| train_02534 | 2noo-A | homolog of CB513 cb513_00420 (2olb-A): identity 0.25, alignment 474, E=4.9e-28 |
| train_02677 | 2om6-B | homolog of CB513 cb513_00193 (1jud-A): identity 0.29, alignment 127, E=8.5e-06 |
| train_02687 | 2opi-B | homolog of CB513 cb513_00127 (1fua-A): identity 0.30, alignment 138, E=1.9e-11 |
| train_02703 | 2ov0-A | homolog of CB513 cb513_00253 (1paz-A): identity 0.37, alignment 79, E=3.3e-05 |
| train_02718 | 2oz8-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.29, alignment 253, E=8.4e-23 |
| train_02766 | 2p75-B | homolog of CB513 cb513_00468 (3pgm-A): identity 0.34, alignment 93, E=3.0e-07 |
| train_02785 | 2pef-A | homolog of CB513 cb513_00509 (9api-B): identity 0.44, alignment 27, E=1.6e-04 |
| train_02793 | 2pgn-B | homolog of CB513 cb513_00270 (1pow-B): identity 0.26, alignment 473, E=2.1e-29 |
| train_02798 | 2pia-A | homolog of CB513 cb513_00131 (1fxi-A): identity 0.32, alignment 62, E=2.9e-06 |
| train_02848 | 2pw8-I | homolog of CB513 cb513_00499 (6hir-A): identity 0.94, alignment 65, E=3.0e-39 |
| train_02885 | 2q4u-A | homolog of CB513 cb513_00455 (3cln-A): identity 0.25, alignment 111, E=2.8e-04 |
| train_02971 | 2qol-A | homolog of CB513 cb513_00189 (1irk-A): identity 0.38, alignment 274, E=5.1e-56 |
| train_02994 | 2qtc-A | homolog of CB513 cb513_00339 (1trk-B): identity 0.28, alignment 160, E=2.2e-05 |
| train_03055 | 2r91-C | homolog of CB513 cb513_00231 (1nal-4): identity 0.25, alignment 280, E=3.6e-17 |
| train_03095 | 2rh7-C | homolog of CB513 cb513_00139 (1gfl-B): identity 0.27, alignment 127, E=1.0e-13 |
| train_03098 | 2rhw-A | homolog of CB513 cb513_00106 (1edd-A): identity 0.33, alignment 119, E=7.3e-09 |
| train_03105 | 2rjo-A | homolog of CB513 cb513_00398 (2gbp-A): identity 0.39, alignment 79, E=5.0e-06 |
| train_03110 | 2rkq-A | homolog of CB513 cb513_00205 (1lba-A): identity 0.26, alignment 111, E=3.8e-06 |
| train_03120 | 2uu8-A | homolog of CB513 cb513_00413 (2ltn-B): identity 0.50, alignment 48, E=1.4e-09 |
| train_03129 | 2uxw-A | homolog of CB513 cb513_00465 (3mdd-B): identity 0.33, alignment 336, E=7.4e-49 |
| train_03291 | 2wg7-A | homolog of CB513 cb513_00476 (4bp2-A): identity 0.50, alignment 35, E=6.4e-04 |
| train_03304 | 2wlr-A | homolog of CB513 cb513_00291 (1rhd-A): identity 0.27, alignment 287, E=2.9e-19 |
| train_03354 | 2x3n-A | homolog of CB513 cb513_00427 (2phh-A): identity 0.26, alignment 318, E=2.4e-10 |
| train_03372 | 2x78-A | homolog of CB513 cb513_00014 (1asw-A): identity 0.34, alignment 92, E=1.5e-04 |
| train_03415 | 2xki-A | homolog of CB513 cb513_00101 (1eca-A): identity 0.33, alignment 85, E=2.1e-06 |
| train_03473 | 2y78-A | homolog of CB513 cb513_00125 (1fkf-A): identity 0.56, alignment 106, E=3.1e-34 |
| train_03593 | 2zcm-A | homolog of CB513 cb513_00444 (2trt-A): identity 0.28, alignment 72, E=3.1e-04 |
| train_03602 | 2zfd-A | homolog of CB513 cb513_00287 (1rec-A): identity 0.30, alignment 115, E=3.7e-08 |
| train_03701 | 3aj7-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.28, alignment 178, E=2.7e-07 |
| train_03715 | 3aml-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.25, alignment 153, E=5.1e-05 |
| train_03747 | 3b0g-A | homolog of CB513 cb513_00399 (2gep-A): identity 0.25, alignment 467, E=5.4e-26 |
| train_03815 | 3bje-B | homolog of CB513 cb513_00105 (1ecp-F): identity 0.28, alignment 120, E=5.1e-06 |
| train_03981 | 3cji-D | homolog of CB513 cb513_00414 (2mev-4): identity 0.40, alignment 65, E=7.0e-06 |
| train_03984 | 3cjm-A | homolog of CB513 cb513_00451 (3blm-A): identity 0.28, alignment 123, E=6.1e-06 |
| train_03996 | 3clm-A | homolog of CB513 cb513_00247 (1onr-B): identity 0.32, alignment 131, E=2.3e-07 |
| train_03999 | 3cm3-A | homolog of CB513 cb513_00352 (1vhr-B): identity 0.29, alignment 112, E=2.2e-07 |
| train_04080 | 3d2l-D | homolog of CB513 cb513_00362 (1xva-B): identity 0.26, alignment 145, E=6.1e-07 |
| train_04145 | 3dhu-D | homolog of CB513 cb513_00051 (1cgu-A): identity 0.26, alignment 296, E=1.9e-13 |
| train_04197 | 3duw-A | homolog of CB513 cb513_00353 (1vid-A): identity 0.27, alignment 165, E=5.9e-13 |
| train_04215 | 3e03-A | homolog of CB513 cb513_00121 (1fdt-A): identity 0.30, alignment 137, E=3.2e-06 |
| train_04232 | 3e3x-A | homolog of CB513 cb513_00079 (1dar-A): identity 0.27, alignment 84, E=9.8e-06 |
| train_04242 | 3e5t-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.29, alignment 129, E=4.5e-13 |
| train_04253 | 3e97-A | homolog of CB513 cb513_00132 (1g6n-A): identity 0.26, alignment 177, E=7.4e-14 |
| train_04271 | 3edf-B | homolog of CB513 cb513_00051 (1cgu-A): identity 0.26, alignment 448, E=8.5e-21 |
| train_04297 | 3ekh-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_04312 | 3eof-A | homolog of CB513 cb513_00241 (1nox-A): identity 0.29, alignment 186, E=5.6e-13 |
| train_04314 | 3ep1-A | homolog of CB513 cb513_00205 (1lba-A): identity 0.34, alignment 92, E=3.2e-07 |
| train_04333 | 3eun-A | homolog of CB513 cb513_00099 (1dur-A): identity 0.40, alignment 57, E=3.0e-09 |
| train_04337 | 3evp-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.96, alignment 143, E=9.0e-90 |
| train_04338 | 3evr-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_04339 | 3evu-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_04361 | 3f1l-B | homolog of CB513 cb513_00121 (1fdt-A): identity 0.28, alignment 167, E=6.2e-09 |
| train_04476 | 3fou-A | homolog of CB513 cb513_00293 (1rie-A): identity 0.36, alignment 75, E=6.1e-06 |
| train_04557 | 3g9m-B | homolog of CB513 cb513_00204 (1lat-B): identity 0.86, alignment 79, E=1.1e-43 |
| train_04563 | 3gb5-A | homolog of CB513 cb513_00241 (1nox-A): identity 0.26, alignment 174, E=4.9e-09 |
| train_04578 | 3gfa-A | homolog of CB513 cb513_00241 (1nox-A): identity 0.39, alignment 57, E=3.2e-07 |
| train_04592 | 3gj0-A | homolog of CB513 cb513_00506 (821p-A): identity 0.26, alignment 161, E=5.9e-14 |
| train_04616 | 3go2-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.25, alignment 354, E=1.1e-11 |
| train_04712 | 3h8g-C | homolog of CB513 cb513_00203 (1lap-A): identity 0.36, alignment 441, E=4.8e-72 |
| train_04721 | 3ha2-A | homolog of CB513 cb513_00283 (1qrd-B): identity 0.32, alignment 73, E=1.2e-06 |
| train_04747 | 3hje-A | homolog of CB513 cb513_00051 (1cgu-A): identity 0.31, alignment 77, E=4.7e-04 |
| train_04767 | 3hno-C | homolog of CB513 cb513_00480 (4pfk-A): identity 0.28, alignment 279, E=1.7e-15 |
| train_04840 | 3i6e-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.27, alignment 268, E=4.5e-18 |
| train_04852 | 3iay-A | homolog of CB513 cb513_00242 (1noz-B): identity 0.28, alignment 107, E=4.3e-06 |
| train_04887 | 3iib-A | homolog of CB513 cb513_00008 (1amp-A): identity 0.34, alignment 110, E=1.1e-09 |
| train_05004 | 3k1z-A | homolog of CB513 cb513_00193 (1jud-A): identity 0.37, alignment 59, E=4.7e-06 |
| train_05014 | 3k59-A | homolog of CB513 cb513_00242 (1noz-B): identity 0.29, alignment 142, E=1.9e-06 |
| train_05074 | 3khe-A | homolog of CB513 cb513_00455 (3cln-A): identity 0.33, alignment 149, E=7.8e-16 |
| train_05141 | 3l0f-A | homolog of CB513 cb513_00063 (1cpc-L): identity 0.26, alignment 117, E=8.4e-11 |
| train_05158 | 3l4n-A | homolog of CB513 cb513_00199 (1kte-A): identity 0.28, alignment 86, E=1.4e-08 |
| train_05203 | 3lgd-A | homolog of CB513 cb513_00004 (1add-A): identity 0.29, alignment 179, E=8.3e-06 |
| train_05211 | 3lhq-B | homolog of CB513 cb513_00444 (2trt-A): identity 0.27, alignment 97, E=4.1e-05 |
| train_05259 | 3luf-B | homolog of CB513 cb513_00453 (3chy-A): identity 0.27, alignment 125, E=1.5e-05 |
| train_05344 | 3mcw-B | homolog of CB513 cb513_00233 (1nba-C): identity 0.26, alignment 148, E=4.8e-05 |
| train_05387 | 3mmy-G | homolog of CB513 cb513_00152 (1gp2-B): identity 0.26, alignment 124, E=1.5e-06 |
| train_05399 | 3mse-B | homolog of CB513 cb513_00455 (3cln-A): identity 0.28, alignment 97, E=1.5e-04 |
| train_05466 | 3na5-A | homolog of CB513 cb513_00470 (3pmg-B): identity 0.27, alignment 380, E=1.5e-18 |
| train_05586 | 3oa5-B | homolog of CB513 cb513_00072 (1ctn-A): identity 0.29, alignment 442, E=5.5e-36 |
| train_05600 | 3odt-B | homolog of CB513 cb513_00152 (1gp2-B): identity 0.25, alignment 250, E=2.1e-08 |
| train_05632 | 3on2-D | homolog of CB513 cb513_00444 (2trt-A): identity 0.37, alignment 51, E=2.4e-04 |
| train_05652 | 3oqp-A | homolog of CB513 cb513_00233 (1nba-C): identity 0.28, alignment 147, E=1.9e-10 |
| train_05677 | 3ow8-A | homolog of CB513 cb513_00152 (1gp2-B): identity 0.27, alignment 196, E=5.8e-24 |
| train_05729 | 3pfg-A | homolog of CB513 cb513_00362 (1xva-B): identity 0.29, alignment 156, E=1.1e-07 |
| train_05758 | 3poh-A | homolog of CB513 cb513_00394 (2ebn-A): identity 0.42, alignment 276, E=3.3e-52 |
| train_05818 | 3qc0-A | homolog of CB513 cb513_00488 (4xia-B): identity 0.29, alignment 149, E=6.8e-04 |
| train_05822 | 3qdl-D | homolog of CB513 cb513_00241 (1nox-A): identity 0.26, alignment 177, E=9.6e-08 |
| train_05824 | 3qex-A | homolog of CB513 cb513_00242 (1noz-B): identity 0.62, alignment 391, E=9.8e-167 |
| train_05830 | 3qft-A | homolog of CB513 cb513_00126 (1fnd-A): identity 0.28, alignment 225, E=4.1e-14 |
| train_05836 | 3qho-C | homolog of CB513 cb513_00102 (1ece-B): identity 0.45, alignment 346, E=7.6e-89 |
| train_05884 | 3qxc-A | homolog of CB513 cb513_00097 (1dts-A): identity 0.29, alignment 156, E=4.7e-09 |
| train_05887 | 3qy9-A | homolog of CB513 cb513_00083 (1dih-A): identity 0.37, alignment 189, E=3.5e-28 |
| train_05906 | 3r41-A | homolog of CB513 cb513_00106 (1edd-A): identity 0.25, alignment 127, E=1.6e-05 |
| train_05919 | 3r7a-B | homolog of CB513 cb513_00468 (3pgm-A): identity 0.30, alignment 92, E=9.2e-08 |
| train_05934 | 3ren-B | homolog of CB513 cb513_00044 (1cem-A): identity 0.26, alignment 184, E=8.1e-10 |
| train_05948 | 3rh2-A | homolog of CB513 cb513_00444 (2trt-A): identity 0.40, alignment 47, E=7.4e-05 |
| train_05979 | 3rot-B | homolog of CB513 cb513_00398 (2gbp-A): identity 0.26, alignment 148, E=5.0e-06 |
| train_06035 | 3s8g-B | homolog of CB513 cb513_00076 (1cyx-A): identity 0.31, alignment 68, E=1.9e-06 |
| train_06057 | 3sg2-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_06058 | 3sg3-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_06059 | 3sg4-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=6.6e-90 |
| train_06060 | 3sg5-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=6.6e-90 |
| train_06061 | 3sg6-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_06062 | 3sg7-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.60, alignment 255, E=9.9e-91 |
| train_06103 | 3sto-A | homolog of CB513 cb513_00508 (9api-A): identity 0.25, alignment 249, E=2.1e-19 |
| train_06107 | 3swo-D | homolog of CB513 cb513_00465 (3mdd-B): identity 0.25, alignment 378, E=7.7e-33 |
| train_06140 | 3t9w-A | homolog of CB513 cb513_00013 (1aoz-B): identity 0.31, alignment 88, E=2.2e-05 |
| train_06170 | 3tjr-A | homolog of CB513 cb513_00121 (1fdt-A): identity 0.29, alignment 178, E=1.1e-12 |
| train_06267 | 3ufb-A | homolog of CB513 cb513_00373 (2adm-B): identity 0.27, alignment 224, E=3.2e-04 |
| train_06270 | 3uff-B | homolog of CB513 cb513_00275 (1ptr-A): identity 0.98, alignment 50, E=4.1e-33 |
| train_06283 | 3ulj-B | homolog of CB513 cb513_00223 (1mjc-A): identity 0.43, alignment 63, E=5.0e-11 |
| train_06337 | 3v7d-D | homolog of CB513 cb513_00152 (1gp2-B): identity 0.25, alignment 222, E=6.9e-16 |
| train_06385 | 3vot-B | homolog of CB513 cb513_00392 (2dln-A): identity 0.29, alignment 225, E=9.4e-11 |
| train_06406 | 3vuq-C | homolog of CB513 cb513_00444 (2trt-A): identity 0.32, alignment 53, E=2.4e-04 |
| train_06477 | 3wlc-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.60, alignment 255, E=2.6e-90 |
| train_06540 | 3zdb-A | homolog of CB513 cb513_00322 (1taq-A): identity 0.33, alignment 243, E=3.2e-28 |
| train_06610 | 4a29-A | homolog of CB513 cb513_00262 (1pii-A): identity 0.32, alignment 197, E=4.8e-22 |
| train_06648 | 4aee-B | homolog of CB513 cb513_00051 (1cgu-A): identity 0.28, alignment 416, E=7.8e-33 |
| train_06841 | 4cf4-B | homolog of CB513 cb513_00193 (1jud-A): identity 0.25, alignment 145, E=1.4e-08 |
| train_06852 | 4ci7-B | homolog of CB513 cb513_00511 (9pap-A): identity 0.25, alignment 196, E=2.8e-06 |
| train_06858 | 4cit-A | homolog of CB513 cb513_00356 (1vnc-A): identity 0.29, alignment 171, E=1.6e-06 |
| train_06861 | 4ck4-B | homolog of CB513 cb513_00111 (1epb-B): identity 0.26, alignment 132, E=1.4e-05 |
| train_06965 | 4dt4-A | homolog of CB513 cb513_00125 (1fkf-A): identity 0.34, alignment 61, E=2.0e-08 |
| train_06970 | 4dwd-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.28, alignment 279, E=5.0e-22 |
| train_07053 | 4eqs-B | homolog of CB513 cb513_00336 (1trb-A): identity 0.26, alignment 221, E=8.0e-11 |
| train_07057 | 4ery-A | homolog of CB513 cb513_00152 (1gp2-B): identity 0.29, alignment 241, E=1.0e-29 |
| train_07065 | 4etn-A | homolog of CB513 cb513_00267 (1pnt-A): identity 0.28, alignment 127, E=1.7e-07 |
| train_07100 | 4f6m-A | homolog of CB513 cb513_00343 (1ubd-C): identity 0.34, alignment 80, E=8.5e-10 |
| train_07135 | 4fin-B | homolog of CB513 cb513_00141 (1gky-A): identity 0.25, alignment 109, E=5.3e-04 |
| train_07281 | 4gy7-A | homolog of CB513 cb513_00198 (1krc-B): identity 0.53, alignment 106, E=9.3e-29 |
| train_07313 | 4h6x-A | homolog of CB513 cb513_00301 (1s01-A): identity 0.27, alignment 183, E=9.2e-05 |
| train_07469 | 4ik1-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.75, alignment 193, E=4.8e-90 |
| train_07470 | 4ik5-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.60, alignment 255, E=9.9e-91 |
| train_07471 | 4ik8-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.60, alignment 255, E=9.9e-91 |
| train_07488 | 4inz-B | homolog of CB513 cb513_00106 (1edd-A): identity 0.32, alignment 116, E=3.0e-12 |
| train_07642 | 4kdr-A | homolog of CB513 cb513_00362 (1xva-B): identity 0.26, alignment 110, E=5.6e-04 |
| train_07646 | 4kem-A | homolog of CB513 cb513_00226 (1mns-A): identity 0.28, alignment 259, E=1.5e-16 |
| train_07698 | 4kyq-A | homolog of CB513 cb513_00352 (1vhr-B): identity 0.26, alignment 127, E=1.7e-09 |
| train_07752 | 4lg8-A | homolog of CB513 cb513_00152 (1gp2-B): identity 0.25, alignment 213, E=4.3e-20 |
| train_07766 | 4lks-A | homolog of CB513 cb513_00147 (1gog-A): identity 0.34, alignment 67, E=1.5e-05 |
| train_07768 | 4lld-B | homolog of CB513 cb513_00505 (7rsa-A): identity 1.00, alignment 20, E=8.7e-05 |
| train_07780 | 4lpq-A | homolog of CB513 cb513_00206 (1lbu-A): identity 0.46, alignment 57, E=5.4e-04 |
| train_07813 | 4lyy-B | homolog of CB513 cb513_00168 (1hmp-B): identity 0.34, alignment 157, E=6.5e-22 |
| train_07873 | 4me9-A | homolog of CB513 cb513_00444 (2trt-A): identity 0.37, alignment 54, E=3.1e-05 |
| train_07916 | 4mt2-A | contains CB513 cb513_00229 (1mrt-A) as a substring |
| train_07926 | 4muv-A | homolog of CB513 cb513_00132 (1g6n-A): identity 0.29, alignment 89, E=2.1e-05 |
| train_07983 | 4ndj-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.94, alignment 238, E=1.8e-152 |
| train_07984 | 4ndk-B | homolog of CB513 cb513_00139 (1gfl-B): identity 0.93, alignment 238, E=4.3e-151 |
| train_08024 | 4npd-A | homolog of CB513 cb513_00119 (1fc2-C): identity 0.88, alignment 58, E=1.5e-27 |
| train_08062 | 4o32-A | homolog of CB513 cb513_00327 (1thx-A): identity 0.26, alignment 85, E=2.7e-11 |
| train_08115 | 4omf-G | homolog of CB513 cb513_00099 (1dur-A): identity 0.34, alignment 47, E=8.2e-04 |
| train_08119 | 4oo4-B | homolog of CB513 cb513_00327 (1thx-A): identity 0.27, alignment 77, E=4.1e-06 |
| train_08211 | 4pf4-A | homolog of CB513 cb513_00189 (1irk-A): identity 0.25, alignment 255, E=1.7e-11 |
| train_08217 | 4phj-B | homolog of CB513 cb513_00455 (3cln-A): identity 0.30, alignment 92, E=6.1e-05 |
| train_08293 | 4q4w-4 | homolog of CB513 cb513_00483 (4rhv-4): identity 0.63, alignment 68, E=8.6e-23 |
| train_08324 | 4qbu-A | homolog of CB513 cb513_00224 (1mla-A): identity 0.25, alignment 291, E=1.3e-09 |
| train_08379 | 4qtt-D | homolog of CB513 cb513_00362 (1xva-B): identity 0.36, alignment 75, E=3.2e-04 |
| train_08422 | 4r81-A | homolog of CB513 cb513_00283 (1qrd-B): identity 0.34, alignment 122, E=2.3e-14 |
| train_08434 | 4rd4-A | homolog of CB513 cb513_00109 (1eft-A): identity 0.32, alignment 263, E=5.6e-27 |
| train_08537 | 4tmx-A | homolog of CB513 cb513_00109 (1eft-A): identity 0.30, alignment 123, E=2.3e-04 |
| train_08573 | 4u1e-I | homolog of CB513 cb513_00152 (1gp2-B): identity 0.26, alignment 177, E=3.0e-11 |
| train_08576 | 4u3y-B | homolog of CB513 cb513_00189 (1irk-A): identity 0.26, alignment 205, E=6.8e-13 |
| train_08630 | 4uho-A | homolog of CB513 cb513_00390 (2dkb-A): identity 0.31, alignment 356, E=2.2e-37 |
| train_08640 | 4uob-A | homolog of CB513 cb513_00371 (2abk-A): identity 0.32, alignment 216, E=2.3e-16 |
| train_08656 | 4uul-B | homolog of CB513 cb513_00386 (2cmd-A): identity 0.27, alignment 205, E=1.2e-06 |
| train_08760 | 4wpk-A | homolog of CB513 cb513_00346 (1udh-A): identity 0.44, alignment 167, E=7.3e-26 |
| train_08768 | 4wtp-A | homolog of CB513 cb513_00140 (1ghs-B): identity 0.28, alignment 158, E=6.1e-04 |
| train_08866 | 4xom-B | homolog of CB513 cb513_00241 (1nox-A): identity 0.38, alignment 66, E=4.8e-05 |
| train_09042 | 4z9h-A | homolog of CB513 cb513_00379 (2asr-A): identity 0.99, alignment 142, E=9.3e-87 |
| train_09068 | 4ze8-D | homolog of CB513 cb513_00420 (2olb-A): identity 0.25, alignment 461, E=4.8e-25 |
| train_09132 | 4zx2-A | homolog of CB513 cb513_00124 (1fjm-B): identity 0.46, alignment 232, E=2.8e-62 |
| train_09136 | 4zyl-B | homolog of CB513 cb513_00453 (3chy-A): identity 0.28, alignment 71, E=7.1e-05 |
| train_09165 | 5a7g-A | homolog of CB513 cb513_00337 (1trh-A): identity 0.33, alignment 319, E=5.2e-35 |
| train_09201 | 5am2-A | homolog of CB513 cb513_00106 (1edd-A): identity 0.33, alignment 131, E=8.9e-15 |
| train_09283 | 5b6d-A | homolog of CB513 cb513_00445 (2tsc-A): identity 0.28, alignment 167, E=3.1e-12 |
| train_09298 | 5bnz-A | homolog of CB513 cb513_00142 (1gln-A): identity 0.26, alignment 294, E=1.9e-10 |
| train_09303 | 5bov-C | homolog of CB513 cb513_00106 (1edd-A): identity 0.27, alignment 121, E=1.2e-05 |
| train_09304 | 5bow-A | homolog of CB513 cb513_00411 (2i1b-A): identity 0.30, alignment 81, E=1.5e-05 |
| train_09315 | 5bql-C | homolog of CB513 cb513_00139 (1gfl-B): identity 0.32, alignment 139, E=1.3e-16 |
| train_09392 | 5cdv-A | homolog of CB513 cb513_00054 (1chm-B): identity 0.27, alignment 271, E=6.7e-22 |
| train_09463 | 5czl-A | homolog of CB513 cb513_00044 (1cem-A): identity 0.31, alignment 180, E=2.7e-12 |
| train_09513 | 5djh-A | homolog of CB513 cb513_00407 (2hhm-B): identity 0.28, alignment 89, E=3.8e-04 |
| train_09550 | 5dze-A | homolog of CB513 cb513_00064 (1cpn-A): identity 0.35, alignment 148, E=3.1e-23 |
| train_09568 | 5e57-A | homolog of CB513 cb513_00444 (2trt-A): identity 0.33, alignment 79, E=9.2e-07 |
| train_09606 | 5eju-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.28, alignment 202, E=3.8e-17 |
| train_09624 | 5eq7-A | homolog of CB513 cb513_00407 (2hhm-B): identity 0.25, alignment 258, E=2.8e-16 |
| train_09642 | 5ew0-A | homolog of CB513 cb513_00366 (1znb-B): identity 0.28, alignment 206, E=5.3e-18 |
| train_09649 | 5exc-O | homolog of CB513 cb513_00139 (1gfl-B): identity 0.28, alignment 140, E=9.7e-11 |
| train_09678 | 5f72-T | homolog of CB513 cb513_00120 (1fdl-H): identity 0.49, alignment 121, E=1.0e-29 |
| train_09690 | 5f9k-B | homolog of CB513 cb513_00098 (1dup-A): identity 0.38, alignment 131, E=2.9e-17 |
| train_09710 | 5fcs-H | homolog of CB513 cb513_00120 (1fdl-H): identity 0.46, alignment 118, E=1.0e-29 |
| train_09711 | 5fcs-L | homolog of CB513 cb513_00120 (1fdl-H): identity 0.50, alignment 125, E=1.5e-32 |
| train_09724 | 5fgu-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.98, alignment 228, E=4.3e-151 |
| train_09735 | 5fky-B | homolog of CB513 cb513_00281 (1qbb-A): identity 0.27, alignment 117, E=1.2e-04 |
| train_09741 | 5fmo-S | homolog of CB513 cb513_00029 (1bmv-1): identity 0.43, alignment 151, E=1.9e-35 |
| train_09751 | 5fpw-A | homolog of CB513 cb513_00511 (9pap-A): identity 0.25, alignment 226, E=2.2e-09 |
| train_09765 | 5fvg-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.31, alignment 139, E=1.4e-20 |
| train_09799 | 5gjh-A | homolog of CB513 cb513_00018 (1aya-B): identity 0.45, alignment 76, E=5.7e-15 |
| train_09804 | 5gkv-A | homolog of CB513 cb513_00363 (1xx2-A): identity 0.26, alignment 160, E=1.5e-04 |
| train_09837 | 5gut-A | homolog of CB513 cb513_00169 (1hmy-A): identity 0.26, alignment 183, E=8.5e-11 |
| train_09846 | 5gy3-A | homolog of CB513 cb513_00044 (1cem-A): identity 0.31, alignment 171, E=1.5e-10 |
| train_09848 | 5gz3-B | homolog of CB513 cb513_00137 (1gd1-O): identity 0.45, alignment 40, E=5.3e-04 |
| train_09876 | 5h88-B | homolog of CB513 cb513_00139 (1gfl-B): identity 0.27, alignment 237, E=3.4e-20 |
| train_09890 | 5hdm-B | homolog of CB513 cb513_00390 (2dkb-A): identity 0.25, alignment 319, E=7.0e-23 |
| train_09907 | 5hlo-B | homolog of CB513 cb513_00462 (3icb-A): identity 0.38, alignment 65, E=3.8e-08 |
| train_09919 | 5hrp-A | homolog of CB513 cb513_00014 (1asw-A): identity 0.35, alignment 105, E=1.3e-05 |
| train_09920 | 5hrr-A | homolog of CB513 cb513_00014 (1asw-A): identity 0.36, alignment 105, E=4.0e-06 |
| train_09921 | 5hrs-A | homolog of CB513 cb513_00014 (1asw-A): identity 0.35, alignment 105, E=9.8e-06 |
| train_09922 | 5hsf-A | homolog of CB513 cb513_00120 (1fdl-H): identity 0.37, alignment 86, E=5.4e-06 |
| train_09934 | 5hwa-A | homolog of CB513 cb513_00053 (1chk-B): identity 0.41, alignment 51, E=1.2e-04 |
| train_09936 | 5hwh-A | homolog of CB513 cb513_00233 (1nba-C): identity 0.27, alignment 169, E=8.6e-06 |
| train_09962 | 5i77-A | homolog of CB513 cb513_00102 (1ece-B): identity 0.25, alignment 216, E=2.5e-04 |
| train_09963 | 5i78-A | homolog of CB513 cb513_00102 (1ece-B): identity 0.25, alignment 216, E=2.5e-04 |
| train_09964 | 5i79-B | homolog of CB513 cb513_00102 (1ece-B): identity 0.25, alignment 216, E=2.5e-04 |
| train_09977 | 5idn-A | homolog of CB513 cb513_00189 (1irk-A): identity 0.26, alignment 209, E=2.1e-06 |
| train_10008 | 5iu4-A | homolog of CB513 cb513_00368 (256b-A): identity 0.98, alignment 105, E=4.1e-59 |
| train_10037 | 5j3n-B | homolog of CB513 cb513_00139 (1gfl-B): identity 0.94, alignment 237, E=6.5e-152 |
| train_10136 | 5k5i-A | homolog of CB513 cb513_00343 (1ubd-C): identity 0.33, alignment 78, E=5.6e-09 |
| train_10174 | 5kia-A | homolog of CB513 cb513_00507 (8adh-A): identity 0.25, alignment 339, E=3.0e-19 |
| train_10184 | 5knl-E | homolog of CB513 cb513_00345 (1ubq-A): identity 0.87, alignment 76, E=5.8e-40 |
| train_10189 | 5krt-A | homolog of CB513 cb513_00014 (1asw-A): identity 0.33, alignment 105, E=1.3e-05 |
| train_10208 | 5l02-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=2.4e-18 |
| train_10210 | 5l05-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.27, alignment 322, E=2.6e-20 |
| train_10279 | 5lnn-A | homolog of CB513 cb513_00453 (3chy-A): identity 0.29, alignment 80, E=2.8e-05 |
| train_10285 | 5lqb-L | homolog of CB513 cb513_00324 (1tcr-A): identity 0.31, alignment 147, E=2.2e-07 |
| train_10313 | 5lxz-B | homolog of CB513 cb513_00147 (1gog-A): identity 0.27, alignment 385, E=1.4e-21 |
| train_10353 | 5mal-A | homolog of CB513 cb513_00112 (1ese-A): identity 0.26, alignment 180, E=1.6e-06 |
| train_10411 | 5mwb-A | homolog of CB513 cb513_00107 (1edm-C): identity 0.48, alignment 35, E=1.4e-10 |
| train_10432 | 5nak-B | homolog of CB513 cb513_00427 (2phh-A): identity 0.25, alignment 333, E=3.7e-08 |
| train_10502 | 5sx1-A | homolog of CB513 cb513_00389 (2cyp-A): identity 0.27, alignment 322, E=2.6e-20 |
| train_10503 | 5sx2-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=9.6e-19 |
| train_10504 | 5sxr-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=9.6e-19 |
| train_10505 | 5sxt-A | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=6.5e-20 |
| train_10506 | 5sxw-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=5.8e-18 |
| train_10507 | 5sxx-B | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=1.6e-19 |
| train_10511 | 5syh-A | homolog of CB513 cb513_00389 (2cyp-A): identity 0.26, alignment 322, E=7.1e-19 |
| train_10512 | 5syi-A | homolog of CB513 cb513_00389 (2cyp-A): identity 0.27, alignment 322, E=2.0e-20 |
| train_10532 | 5t5i-F | homolog of CB513 cb513_00099 (1dur-A): identity 0.35, alignment 57, E=5.4e-08 |
| train_10535 | 5t5l-A | homolog of CB513 cb513_00413 (2ltn-B): identity 0.50, alignment 42, E=1.6e-06 |
| train_10551 | 5tde-A | homolog of CB513 cb513_00303 (1scu-E): identity 0.25, alignment 260, E=1.6e-04 |
| train_10553 | 5tee-A | homolog of CB513 cb513_00152 (1gp2-B): identity 0.25, alignment 167, E=8.7e-08 |
| train_10588 | 5tv7-B | homolog of CB513 cb513_00206 (1lbu-A): identity 0.42, alignment 69, E=1.0e-08 |
| train_10611 | 5u4s-B | homolog of CB513 cb513_00121 (1fdt-A): identity 0.28, alignment 193, E=1.1e-10 |
| train_10645 | 5uiw-A | homolog of CB513 cb513_00484 (4rxn-A): identity 0.94, alignment 54, E=1.1e-34 |
| train_10654 | 5um7-A | homolog of CB513 cb513_00327 (1thx-A): identity 0.38, alignment 40, E=4.4e-04 |
| train_10693 | 5v5h-B | homolog of CB513 cb513_00345 (1ubq-A): identity 0.91, alignment 75, E=2.8e-39 |
| train_10694 | 5v5i-B | homolog of CB513 cb513_00345 (1ubq-A): identity 0.91, alignment 76, E=7.3e-39 |
| train_10709 | 5vfc-A | homolog of CB513 cb513_00152 (1gp2-B): identity 0.29, alignment 241, E=1.0e-29 |
| train_10722 | 5vjy-E | homolog of CB513 cb513_00098 (1dup-A): identity 0.39, alignment 131, E=3.2e-18 |
| train_10739 | 5w1u-B | homolog of CB513 cb513_00337 (1trh-A): identity 0.38, alignment 222, E=2.4e-29 |
| train_10836 | 5y8r-A | homolog of CB513 cb513_00139 (1gfl-B): identity 0.27, alignment 157, E=1.3e-10 |
