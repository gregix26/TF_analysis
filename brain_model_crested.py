# Checking CRESted pre-trained model

import crested

# Download DeepHumanCortex1 — this pulls weights from the CREsted model repo
model_path, output_names = crested.get_model("DeepHumanCortex1")

# Load the model
model = crested.utils.load_model(model_path)

# Check what cell types it covers
print(output_names)
print(model.input_shape)