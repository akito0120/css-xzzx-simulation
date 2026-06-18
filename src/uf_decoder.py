import pathlib
import numpy as np
import scipy.sparse as sp
import stim
from sinter import Decoder
from beliefmatching import detector_error_model_to_check_matrices
from ldpc.union_find_decoder import UnionFindDecoder

class UnionFindSinterDecoder(Decoder):
    def __init__(self, uf_method: str = "peeling"):
        self.uf_method = uf_method

    def decode_via_files(
        self,
        *,
        num_shots: int,
        num_dets: int,
        num_obs: int,
        dem_path: pathlib.Path,
        dets_b8_in_path: pathlib.Path,
        obs_predictions_b8_out_path: pathlib.Path,
        tmp_dir: pathlib.Path,
    ) -> None:
        dem = stim.DetectorErrorModel.from_file(dem_path)
        matrices = detector_error_model_to_check_matrices(dem)

        # Edge-level matrices have columns of weight <= 2, which Union-Find needs
        check_matrix = matrices.edge_check_matrix
        observables_matrix = matrices.edge_observables_matrix

        # Per-edge prior: largest prior among the hyperedges decomposing into it
        priors = np.asarray(matrices.priors, dtype=float)
        edge_prior = np.asarray(
            (matrices.hyperedge_to_edge_matrix @ sp.diags(priors)).max(axis=1).todense()
        ).ravel()
        edge_prior = np.clip(edge_prior, 1e-12, 0.5)
        llrs = np.log((1.0 - edge_prior) / edge_prior)

        uf = UnionFindDecoder(check_matrix, uf_method=self.uf_method)

        shots = stim.read_shot_data_file(
            path=dets_b8_in_path,
            format="b8",
            num_detectors=dem.num_detectors,
            bit_packed=False,
        )

        predictions = np.zeros((num_shots, num_obs), dtype=np.uint8)
        for i, syndrome in enumerate(shots):
            correction = uf.decode(syndrome.astype(np.uint8), llrs)
            predictions[i] = (observables_matrix @ correction) % 2

        stim.write_shot_data_file(
            data=predictions,
            path=obs_predictions_b8_out_path,
            format="b8",
            num_observables=dem.num_observables,
        )
