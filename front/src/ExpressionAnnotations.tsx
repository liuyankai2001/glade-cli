import { useId } from "react";
import { palette } from "./palette";
import type { Feature } from "./types";
import {
  annotationAngle, annotationArc, annotationColors, annotationName, annotationNames,
  annotationPoint, annotationRole, annotationRoles, annotationTooltip, annotationTracks,
  geneLabels, genePath, textArc, textWidth,
} from "./annotationLayout";
import type { AnnotationModel, ReservedLabel } from "./annotationLayout";

export function ExpressionAnnotations({ model, length, reserved }: {
  model: AnnotationModel; length: number; reserved: ReservedLabel[];
}) {
  const prefix = useId();
  const labels = geneLabels(model.items, length, reserved);
  return <g className="expression-annotations" onMouseDown={event => event.stopPropagation()} onClick={event => event.stopPropagation()}>
    <defs>
      {labels.filter(label => label.fits).map(({ item, start, end }) =>
        <path key={item.key} id={`${prefix}-${item.key}`} d={textArc(103, start, end)} />)}
      {model.groups.map(group => <path key={group.key} id={`${prefix}-cassette-${group.key}`}
        d={textArc(120, annotationAngle(group.feature.start_bp, length), annotationAngle(group.feature.end_bp + 1, length))} />)}
    </defs>
    {model.groups.map((group, index) => {
      const start = annotationAngle(group.feature.start_bp, length), end = annotationAngle(group.feature.end_bp + 1, length);
      return <g key={group.key} data-cassette-group={group.key}>
        <title>{group.name} · {group.feature.label} · {group.feature.start_bp}–{group.feature.end_bp} bp</title>
        <path className="cassette-band" d={annotationArc(122, start, end)} fill="none"
          stroke={index % 2 ? palette.expression : palette.gene} strokeWidth="2" />
        {(end - start) * 120 > textWidth(group.name) + 10 && <text className="cassette-name">
          <textPath href={`#${prefix}-cassette-${group.key}`} startOffset="50%" textAnchor="middle">{group.name}</textPath>
        </text>}
      </g>;
    })}
    {model.items.map(item => {
      const { feature, role } = item;
      const radius = annotationTracks[role];
      const mid = (annotationAngle(feature.start_bp, length) + annotationAngle(feature.end_bp + 1, length)) / 2;
      const [x, y] = annotationPoint(radius, mid);
      const directed = feature.strand === 1 || feature.strand === -1;
      return <g key={item.key} data-annotation-role={role} data-start-bp={feature.start_bp}
        data-end-bp={feature.end_bp} data-track-radius={radius} className="annotation-item">
        <title>{annotationTooltip(feature)}</title>
        {role === "gene" ? <path data-testid={directed ? "strand-arrow" : undefined}
          data-strand={feature.strand} data-start-bp={feature.start_bp} data-end-bp={feature.end_bp}
          d={genePath(feature, length)} fill={annotationColors.gene} stroke="#0d1117" strokeWidth="0.6" />
          : <g className="annotation-marker-anchor" transform={`translate(${x} ${y})`}>
            <circle r="5" fill="#0d1117" />
            <text className="annotation-symbol" textAnchor="middle" y="2.6" style={{ fill: annotationColors[role] }}>
              {role === "promoter" ? "P" : role === "rbs" ? "R" : "T"}
            </text>
          </g>}
      </g>;
    })}
    {labels.map(({ item, fits, mid, leader }) => {
      if (fits) return <text key={item.key} className="annotation-gene-name" data-gene-label={item.key} data-label-placement="arc" dy="2.8">
        <title>{annotationTooltip(item.feature)}</title>
        <textPath href={`#${prefix}-${item.key}`} startOffset="50%" textAnchor="middle">{item.name}</textPath>
      </text>;
      if (!leader) return null;
      const [x, y] = annotationPoint(110, mid);
      return <g key={item.key} className="annotation-gene-leader" pointerEvents="none">
        <path d={`M ${x} ${y} L ${leader.x - leader.side * 40} ${leader.y} L ${leader.x - leader.side * 4} ${leader.y}`}
          fill="none" stroke={annotationColors.gene} strokeWidth="1" />
        <text data-gene-label={item.key} data-label-placement="leader" className="annotation-leader-name"
          x={leader.x} y={leader.y - 4} textAnchor={leader.side > 0 ? "end" : "start"}>
          <title>{annotationTooltip(item.feature)}</title>{leader.name}
        </text>
      </g>;
    })}
  </g>;
}

function FeatureRow({ feature }: { feature: Feature }) {
  const role = annotationRole(feature);
  const name = annotationName(feature);
  return <li title={annotationTooltip(feature)}>
    <i style={{ background: role ? annotationColors[role] : palette[feature.kind === "gene" ? "linker" : feature.kind] || "#77839a" }} />{" "}
    <span>{name}</span> · {feature.start_bp}–{feature.end_bp} bp{" "}
    {role && (feature.strand === 1 ? "→" : feature.strand === -1 ? "←" : "方向未标注")}
    {name !== feature.label && <small className="annotation-original-name">{feature.label}</small>}
  </li>;
}

export function AnnotationDetails({ model, features }: { model: AnnotationModel; features: Feature[] }) {
  return <>
    <div className="annotation-key" data-testid="annotation-key" aria-label="内部元件图例">
      {annotationRoles.map(role => <span key={role}>
        <i style={{ background: annotationColors[role] }} />{annotationNames[role]}
      </span>)}
    </div>
    <details className="annotation-details">
      <summary>注释详情 · {features.length} 项</summary>
      <div className="annotation-scroll">
        {model.groups.map(group => <section key={group.key} data-cassette-details={group.key}>
          <h3>{group.name}</h3>
          <small>{group.feature.label} · {group.feature.start_bp}–{group.feature.end_bp} bp</small>
          <ul className="feature-legend">{group.items.map(item => <FeatureRow key={item.key} feature={item.feature} />)}</ul>
        </section>)}
        {model.remaining.length > 0 && <>
          {model.groups.length > 0 && <h3>其他注释</h3>}
          <ul className="feature-legend">{model.remaining.map((feature, index) => <FeatureRow key={index} feature={feature} />)}</ul>
        </>}
      </div>
    </details>
  </>;
}
