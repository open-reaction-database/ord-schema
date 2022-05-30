# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Text generation script for Reaction messages.

This script is meant to convert a Reaction message into a full-fledged text
paragraph that would be appropriate for inclusion in a Supplemental Information
document for a publication.

Example usage:
* For normal operation from a pbtxt
  $ python generate_text.py --input_reaction=reaction.pb --input_format=pbtxt
      --type text
"""
# pylint: disable=import-error
from absl import app
from absl import flags

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text

FLAGS = flags.FLAGS
flags.DEFINE_string("input", None, "File containing a Reaction message.")
flags.DEFINE_enum("output_type", "text", ["text", "html"], "Text or HTML output format.")
flags.DEFINE_string("output", None, "Filename for output Dataset.")


def main(argv):
    del argv  # Only used by app.run().
    reaction = message_helpers.load_message(FLAGS.input, reaction_pb2.Reaction)
    if FLAGS.output_type == "html":
        text = generate_text.generate_html(reaction)
    elif FLAGS.output_type == "text":
        text = generate_text.generate_text(reaction)
    else:
        raise ValueError(f"unsupported output_type: {FLAGS.output_type}")
    if FLAGS.output:
        with open(FLAGS.output, "w") as f:
            f.write(text)
    else:
        print(text)


if __name__ == "__main__":
    flags.mark_flag_as_required("input")
    app.run(main)
