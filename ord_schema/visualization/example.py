# Copyright 2021 Open Reaction Database Project Authors
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
"""Example reaction for tuning templates."""

import os

from google.protobuf import text_format

from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text

if __name__ == '__main__':
    reaction = reaction_pb2.Reaction()
    with open(
            os.path.join(os.path.dirname(__file__), 'testdata',
                         'reaction.pbtxt')) as f:
        text_format.Parse(f.read(), reaction)
    with open('example.html', 'w') as f:
        f.write(generate_text.generate_summary(reaction, dataset_id='test'))
