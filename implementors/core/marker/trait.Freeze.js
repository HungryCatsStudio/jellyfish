(function() {var implementors = {};
implementors["jf_plonk"] = [{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/basic/struct.PlonkCircuit.html\" title=\"struct jf_plonk::circuit::basic::PlonkCircuit\">PlonkCircuit</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::basic::PlonkCircuit"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ecc/struct.Point.html\" title=\"struct jf_plonk::circuit::customized::ecc::Point\">Point</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::customized::ecc::Point"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ecc/struct.PointVariable.html\" title=\"struct jf_plonk::circuit::customized::ecc::PointVariable\">PointVariable</a>","synthetic":true,"types":["jf_plonk::circuit::customized::ecc::PointVariable"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/rescue/struct.RescueStateVar.html\" title=\"struct jf_plonk::circuit::customized::rescue::RescueStateVar\">RescueStateVar</a>","synthetic":true,"types":["jf_plonk::circuit::customized::rescue::native::RescueStateVar"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/rescue/struct.RescueNonNativeStateVar.html\" title=\"struct jf_plonk::circuit::customized::rescue::RescueNonNativeStateVar\">RescueNonNativeStateVar</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::customized::rescue::non_native::RescueNonNativeStateVar"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/transcript/struct.RescueTranscriptVar.html\" title=\"struct jf_plonk::circuit::customized::transcript::RescueTranscriptVar\">RescueTranscriptVar</a>&lt;F&gt;","synthetic":true,"types":["jf_plonk::circuit::customized::transcript::RescueTranscriptVar"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ultraplonk/mod_arith/struct.FpElem.html\" title=\"struct jf_plonk::circuit::customized::ultraplonk::mod_arith::FpElem\">FpElem</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::customized::ultraplonk::mod_arith::FpElem"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ultraplonk/mod_arith/struct.FpElemVar.html\" title=\"struct jf_plonk::circuit::customized::ultraplonk::mod_arith::FpElemVar\">FpElemVar</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::customized::ultraplonk::mod_arith::FpElemVar"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ultraplonk/plonk_verifier/struct.BatchProofVar.html\" title=\"struct jf_plonk::circuit::customized::ultraplonk::plonk_verifier::BatchProofVar\">BatchProofVar</a>&lt;F&gt;","synthetic":true,"types":["jf_plonk::circuit::customized::ultraplonk::plonk_verifier::structs::BatchProofVar"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/customized/ultraplonk/plonk_verifier/struct.VerifyingKeyVar.html\" title=\"struct jf_plonk::circuit::customized::ultraplonk::plonk_verifier::VerifyingKeyVar\">VerifyingKeyVar</a>&lt;E&gt;","synthetic":true,"types":["jf_plonk::circuit::customized::ultraplonk::plonk_verifier::VerifyingKeyVar"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.PaddingGate.html\" title=\"struct jf_plonk::circuit::gates::PaddingGate\">PaddingGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::PaddingGate"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.ConstantGate.html\" title=\"struct jf_plonk::circuit::gates::ConstantGate\">ConstantGate</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::gates::ConstantGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.AdditionGate.html\" title=\"struct jf_plonk::circuit::gates::AdditionGate\">AdditionGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::AdditionGate"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.ConstantAdditionGate.html\" title=\"struct jf_plonk::circuit::gates::ConstantAdditionGate\">ConstantAdditionGate</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::gates::ConstantAdditionGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.SubtractionGate.html\" title=\"struct jf_plonk::circuit::gates::SubtractionGate\">SubtractionGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::SubtractionGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.MultiplicationGate.html\" title=\"struct jf_plonk::circuit::gates::MultiplicationGate\">MultiplicationGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::MultiplicationGate"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.ConstantMultiplicationGate.html\" title=\"struct jf_plonk::circuit::gates::ConstantMultiplicationGate\">ConstantMultiplicationGate</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::circuit::gates::ConstantMultiplicationGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.BoolGate.html\" title=\"struct jf_plonk::circuit::gates::BoolGate\">BoolGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::BoolGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.EqualityGate.html\" title=\"struct jf_plonk::circuit::gates::EqualityGate\">EqualityGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::EqualityGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.IoGate.html\" title=\"struct jf_plonk::circuit::gates::IoGate\">IoGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::IoGate"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/circuit/gates/struct.FifthRootGate.html\" title=\"struct jf_plonk::circuit::gates::FifthRootGate\">FifthRootGate</a>","synthetic":true,"types":["jf_plonk::circuit::gates::FifthRootGate"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_plonk/errors/enum.PlonkError.html\" title=\"enum jf_plonk::errors::PlonkError\">PlonkError</a>","synthetic":true,"types":["jf_plonk::errors::PlonkError"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_plonk/errors/enum.SnarkError.html\" title=\"enum jf_plonk::errors::SnarkError\">SnarkError</a>","synthetic":true,"types":["jf_plonk::errors::SnarkError"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_plonk/errors/enum.CircuitError.html\" title=\"enum jf_plonk::errors::CircuitError\">CircuitError</a>","synthetic":true,"types":["jf_plonk::errors::CircuitError"]},{"text":"impl&lt;'a, E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/batch_arg/struct.BatchArgument.html\" title=\"struct jf_plonk::proof_system::batch_arg::BatchArgument\">BatchArgument</a>&lt;'a, E&gt;","synthetic":true,"types":["jf_plonk::proof_system::batch_arg::BatchArgument"]},{"text":"impl&lt;'a, E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/batch_arg/struct.Instance.html\" title=\"struct jf_plonk::proof_system::batch_arg::Instance\">Instance</a>&lt;'a, E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::Fr: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Prepared: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::batch_arg::Instance"]},{"text":"impl&lt;'a, E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/struct.PlonkKzgSnark.html\" title=\"struct jf_plonk::proof_system::PlonkKzgSnark\">PlonkKzgSnark</a>&lt;'a, E&gt;","synthetic":true,"types":["jf_plonk::proof_system::snark::PlonkKzgSnark"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.UniversalSrs.html\" title=\"struct jf_plonk::proof_system::structs::UniversalSrs\">UniversalSrs</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Prepared: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::UniversalSrs"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.Proof.html\" title=\"struct jf_plonk::proof_system::structs::Proof\">Proof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::Fr: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::Proof"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.PlookupProof.html\" title=\"struct jf_plonk::proof_system::structs::PlookupProof\">PlookupProof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::Fr: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::PlookupProof"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.BatchProof.html\" title=\"struct jf_plonk::proof_system::structs::BatchProof\">BatchProof</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::BatchProof"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.ProofEvaluations.html\" title=\"struct jf_plonk::proof_system::structs::ProofEvaluations\">ProofEvaluations</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::ProofEvaluations"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.PlookupEvaluations.html\" title=\"struct jf_plonk::proof_system::structs::PlookupEvaluations\">PlookupEvaluations</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::PlookupEvaluations"]},{"text":"impl&lt;'a, E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.ProvingKey.html\" title=\"struct jf_plonk::proof_system::structs::ProvingKey\">ProvingKey</a>&lt;'a, E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Prepared: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::ProvingKey"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.PlookupProvingKey.html\" title=\"struct jf_plonk::proof_system::structs::PlookupProvingKey\">PlookupProvingKey</a>&lt;E&gt;","synthetic":true,"types":["jf_plonk::proof_system::structs::PlookupProvingKey"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.VerifyingKey.html\" title=\"struct jf_plonk::proof_system::structs::VerifyingKey\">VerifyingKey</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Affine: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G2Prepared: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::VerifyingKey"]},{"text":"impl&lt;E&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/proof_system/structs/struct.PlookupVerifyingKey.html\" title=\"struct jf_plonk::proof_system::structs::PlookupVerifyingKey\">PlookupVerifyingKey</a>&lt;E&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;E as PairingEngine&gt;::G1Affine: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::proof_system::structs::PlookupVerifyingKey"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_plonk/transcript/struct.RescueTranscript.html\" title=\"struct jf_plonk::transcript::RescueTranscript\">RescueTranscript</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_plonk::transcript::rescue::RescueTranscript"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/transcript/struct.SolidityTranscript.html\" title=\"struct jf_plonk::transcript::SolidityTranscript\">SolidityTranscript</a>","synthetic":true,"types":["jf_plonk::transcript::solidity::SolidityTranscript"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_plonk/transcript/struct.StandardTranscript.html\" title=\"struct jf_plonk::transcript::StandardTranscript\">StandardTranscript</a>","synthetic":true,"types":["jf_plonk::transcript::standard::StandardTranscript"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_plonk/enum.PlonkType.html\" title=\"enum jf_plonk::PlonkType\">PlonkType</a>","synthetic":true,"types":["jf_plonk::PlonkType"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_plonk/enum.MergeableCircuitType.html\" title=\"enum jf_plonk::MergeableCircuitType\">MergeableCircuitType</a>","synthetic":true,"types":["jf_plonk::MergeableCircuitType"]}];
implementors["jf_primitives"] = [{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/aead/struct.EncKey.html\" title=\"struct jf_primitives::aead::EncKey\">EncKey</a>","synthetic":true,"types":["jf_primitives::aead::EncKey"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/aead/struct.KeyPair.html\" title=\"struct jf_primitives::aead::KeyPair\">KeyPair</a>","synthetic":true,"types":["jf_primitives::aead::KeyPair"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/aead/struct.Ciphertext.html\" title=\"struct jf_primitives::aead::Ciphertext\">Ciphertext</a>","synthetic":true,"types":["jf_primitives::aead::Ciphertext"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/elgamal/struct.EncKeyVars.html\" title=\"struct jf_primitives::circuit::elgamal::EncKeyVars\">EncKeyVars</a>","synthetic":true,"types":["jf_primitives::circuit::elgamal::EncKeyVars"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/elgamal/struct.ElGamalHybridCtxtVars.html\" title=\"struct jf_primitives::circuit::elgamal::ElGamalHybridCtxtVars\">ElGamalHybridCtxtVars</a>","synthetic":true,"types":["jf_primitives::circuit::elgamal::ElGamalHybridCtxtVars"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/merkle_tree/struct.MerkleNodeVars.html\" title=\"struct jf_primitives::circuit::merkle_tree::MerkleNodeVars\">MerkleNodeVars</a>","synthetic":true,"types":["jf_primitives::circuit::merkle_tree::MerkleNodeVars"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/merkle_tree/struct.MerklePathVars.html\" title=\"struct jf_primitives::circuit::merkle_tree::MerklePathVars\">MerklePathVars</a>","synthetic":true,"types":["jf_primitives::circuit::merkle_tree::MerklePathVars"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/merkle_tree/struct.AccElemVars.html\" title=\"struct jf_primitives::circuit::merkle_tree::AccElemVars\">AccElemVars</a>","synthetic":true,"types":["jf_primitives::circuit::merkle_tree::AccElemVars"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/merkle_tree/struct.AccMemberWitnessVar.html\" title=\"struct jf_primitives::circuit::merkle_tree::AccMemberWitnessVar\">AccMemberWitnessVar</a>","synthetic":true,"types":["jf_primitives::circuit::merkle_tree::AccMemberWitnessVar"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/signature/schnorr/struct.VerKeyVar.html\" title=\"struct jf_primitives::circuit::signature::schnorr::VerKeyVar\">VerKeyVar</a>","synthetic":true,"types":["jf_primitives::circuit::signature::schnorr::VerKeyVar"]},{"text":"impl Freeze for <a class=\"struct\" href=\"jf_primitives/circuit/signature/schnorr/struct.SignatureVar.html\" title=\"struct jf_primitives::circuit::signature::schnorr::SignatureVar\">SignatureVar</a>","synthetic":true,"types":["jf_primitives::circuit::signature::schnorr::SignatureVar"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/commitment/struct.Commitment.html\" title=\"struct jf_primitives::commitment::Commitment\">Commitment</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::commitment::Commitment"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/elgamal/struct.EncKey.html\" title=\"struct jf_primitives::elgamal::EncKey\">EncKey</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::elgamal::EncKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/elgamal/struct.KeyPair.html\" title=\"struct jf_primitives::elgamal::KeyPair\">KeyPair</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::ScalarField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::elgamal::KeyPair"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/elgamal/struct.Ciphertext.html\" title=\"struct jf_primitives::elgamal::Ciphertext\">Ciphertext</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::elgamal::Ciphertext"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_primitives/errors/enum.PrimitivesError.html\" title=\"enum jf_primitives::errors::PrimitivesError\">PrimitivesError</a>","synthetic":true,"types":["jf_primitives::errors::PrimitivesError"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_primitives/merkle_tree/enum.NodePos.html\" title=\"enum jf_primitives::merkle_tree::NodePos\">NodePos</a>","synthetic":true,"types":["jf_primitives::merkle_tree::NodePos"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerklePathNode.html\" title=\"struct jf_primitives::merkle_tree::MerklePathNode\">MerklePathNode</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerklePathNode"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerklePath.html\" title=\"struct jf_primitives::merkle_tree::MerklePath\">MerklePath</a>&lt;F&gt;","synthetic":true,"types":["jf_primitives::merkle_tree::MerklePath"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.NodeValue.html\" title=\"struct jf_primitives::merkle_tree::NodeValue\">NodeValue</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::NodeValue"]},{"text":"impl&lt;F, P&gt; Freeze for <a class=\"enum\" href=\"jf_primitives/merkle_tree/enum.LookupResult.html\" title=\"enum jf_primitives::merkle_tree::LookupResult\">LookupResult</a>&lt;F, P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;P: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::LookupResult"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerkleCommitment.html\" title=\"struct jf_primitives::merkle_tree::MerkleCommitment\">MerkleCommitment</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerkleCommitment"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerkleLeaf.html\" title=\"struct jf_primitives::merkle_tree::MerkleLeaf\">MerkleLeaf</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerkleLeaf"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerkleLeafProof.html\" title=\"struct jf_primitives::merkle_tree::MerkleLeafProof\">MerkleLeafProof</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerkleLeafProof"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"enum\" href=\"jf_primitives/merkle_tree/enum.MerkleFrontier.html\" title=\"enum jf_primitives::merkle_tree::MerkleFrontier\">MerkleFrontier</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerkleFrontier"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.MerkleTree.html\" title=\"struct jf_primitives::merkle_tree::MerkleTree\">MerkleTree</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::MerkleTree"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.FilledMTBuilder.html\" title=\"struct jf_primitives::merkle_tree::FilledMTBuilder\">FilledMTBuilder</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::FilledMTBuilder"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/merkle_tree/struct.AccMemberWitness.html\" title=\"struct jf_primitives::merkle_tree::AccMemberWitness\">AccMemberWitness</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::merkle_tree::AccMemberWitness"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/prf/struct.PRF.html\" title=\"struct jf_primitives::prf::PRF\">PRF</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::prf::PRF"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/prf/struct.PrfKey.html\" title=\"struct jf_primitives::prf::PrfKey\">PrfKey</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::prf::PrfKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/bls/struct.BLSSignatureScheme.html\" title=\"struct jf_primitives::signatures::bls::BLSSignatureScheme\">BLSSignatureScheme</a>&lt;P&gt;","synthetic":true,"types":["jf_primitives::signatures::bls::BLSSignatureScheme"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/bls/struct.BLSVerKey.html\" title=\"struct jf_primitives::signatures::bls::BLSVerKey\">BLSVerKey</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as Bls12Parameters&gt;::Fp: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::bls::BLSVerKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/bls/struct.BLSSignKey.html\" title=\"struct jf_primitives::signatures::bls::BLSSignKey\">BLSSignKey</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;&lt;P as Bls12Parameters&gt;::G1Parameters as ModelParameters&gt;::ScalarField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::bls::BLSSignKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/bls/struct.BLSSignature.html\" title=\"struct jf_primitives::signatures::bls::BLSSignature\">BLSSignature</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as Bls12Parameters&gt;::Fp: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::bls::BLSSignature"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/schnorr/struct.SchnorrSignatureScheme.html\" title=\"struct jf_primitives::signatures::schnorr::SchnorrSignatureScheme\">SchnorrSignatureScheme</a>&lt;P&gt;","synthetic":true,"types":["jf_primitives::signatures::schnorr::SchnorrSignatureScheme"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/schnorr/struct.SignKey.html\" title=\"struct jf_primitives::signatures::schnorr::SignKey\">SignKey</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::schnorr::SignKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/schnorr/struct.VerKey.html\" title=\"struct jf_primitives::signatures::schnorr::VerKey\">VerKey</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::schnorr::VerKey"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/schnorr/struct.KeyPair.html\" title=\"struct jf_primitives::signatures::schnorr::KeyPair\">KeyPair</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::ScalarField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::schnorr::KeyPair"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/signatures/schnorr/struct.Signature.html\" title=\"struct jf_primitives::signatures::schnorr::Signature\">Signature</a>&lt;P&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::BaseField: Freeze,<br>&nbsp;&nbsp;&nbsp;&nbsp;&lt;P as ModelParameters&gt;::ScalarField: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_primitives::signatures::schnorr::Signature"]},{"text":"impl&lt;P&gt; Freeze for <a class=\"struct\" href=\"jf_primitives/vrf/blsvrf/struct.BLSVRFScheme.html\" title=\"struct jf_primitives::vrf::blsvrf::BLSVRFScheme\">BLSVRFScheme</a>&lt;P&gt;","synthetic":true,"types":["jf_primitives::vrf::blsvrf::BLSVRFScheme"]}];
implementors["jf_rescue"] = [{"text":"impl Freeze for <a class=\"enum\" href=\"jf_rescue/errors/enum.RescueError.html\" title=\"enum jf_rescue::errors::RescueError\">RescueError</a>","synthetic":true,"types":["jf_rescue::errors::RescueError"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_rescue/struct.RescueVector.html\" title=\"struct jf_rescue::RescueVector\">RescueVector</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_rescue::RescueVector"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_rescue/struct.RescueMatrix.html\" title=\"struct jf_rescue::RescueMatrix\">RescueMatrix</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_rescue::RescueMatrix"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_rescue/struct.PRP.html\" title=\"struct jf_rescue::PRP\">PRP</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_rescue::PRP"]},{"text":"impl&lt;F&gt; Freeze for <a class=\"struct\" href=\"jf_rescue/struct.Permutation.html\" title=\"struct jf_rescue::Permutation\">Permutation</a>&lt;F&gt; <span class=\"where fmt-newline\">where<br>&nbsp;&nbsp;&nbsp;&nbsp;F: Freeze,&nbsp;</span>","synthetic":true,"types":["jf_rescue::Permutation"]}];
implementors["jf_utils"] = [{"text":"impl Freeze for <a class=\"struct\" href=\"jf_utils/struct.CanonicalBytes.html\" title=\"struct jf_utils::CanonicalBytes\">CanonicalBytes</a>","synthetic":true,"types":["jf_utils::serialize::CanonicalBytes"]},{"text":"impl&lt;T&gt; Freeze for <a class=\"struct\" href=\"jf_utils/struct.TaggedBlob.html\" title=\"struct jf_utils::TaggedBlob\">TaggedBlob</a>&lt;T&gt;","synthetic":true,"types":["jf_utils::serialize::TaggedBlob"]},{"text":"impl Freeze for <a class=\"enum\" href=\"jf_utils/enum.TaggedBlobError.html\" title=\"enum jf_utils::TaggedBlobError\">TaggedBlobError</a>","synthetic":true,"types":["jf_utils::serialize::TaggedBlobError"]}];
if (window.register_implementors) {window.register_implementors(implementors);} else {window.pending_implementors = implementors;}})()