/**
 * Copyright 2026 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// A map key named __proto__ is ordinary protobuf data: decoding it has to produce an entry of
// the map, not a new prototype for it.
import assert from "node:assert/strict";
import { ord } from "../index.js";

const original = ord.Reaction.fromObject({
    reactionId: "r1",
    inputs: { ["__proto__"]: { additionOrder: 7 } },
});
const decoded = ord.Reaction.decode(ord.Reaction.encode(original).finish());

assert.equal(Object.getPrototypeOf(decoded.inputs), Object.prototype);
assert.ok(Object.hasOwn(decoded.inputs, "__proto__"));
assert.equal(decoded.inputs["__proto__"].additionOrder, 7);
assert.equal(ord.Reaction.toObject(decoded).inputs["__proto__"].additionOrder, 7);
console.log("map keys ok");
