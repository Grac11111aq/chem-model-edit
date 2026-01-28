import { cn } from '@/lib/utils'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'

type AtomTableRow = {
  index: number
  symbol: string
  x: number
  y: number
  z: number
  fixed?: boolean
}

type AtomTableProps = {
  rows: Array<AtomTableRow>
  selectedIndices?: Set<number>
  fixedIndices?: Set<number>
  onRowClick?: (index: number) => void
  selectionEnabled?: boolean
  digits?: number
  emptyText?: string
  containerClassName?: string
  stickyHeader?: boolean
  showIndex?: boolean
}

const formatNumber = (value: number, digits: number) => value.toFixed(digits)

export function AtomTable({
  rows,
  selectedIndices,
  fixedIndices,
  onRowClick,
  selectionEnabled = true,
  digits = 3,
  emptyText = 'No atoms.',
  containerClassName,
  stickyHeader = false,
  showIndex = true,
}: AtomTableProps) {
  const canSelect = Boolean(onRowClick) && selectionEnabled
  const colSpan = showIndex ? 5 : 4

  return (
    <Table
      className="text-xs"
      containerClassName={cn('h-full', containerClassName)}
    >
      <TableHeader className={cn(stickyHeader && 'sticky top-0 bg-slate-50')}>
        <TableRow>
          {showIndex ? <TableHead>#</TableHead> : null}
          <TableHead>El</TableHead>
          <TableHead>X</TableHead>
          <TableHead>Y</TableHead>
          <TableHead>Z</TableHead>
        </TableRow>
      </TableHeader>
      <TableBody>
        {rows.length ? (
          rows.map((row) => {
            const isFixed = fixedIndices?.has(row.index) ?? row.fixed ?? false
            const isSelected = selectedIndices?.has(row.index) ?? false
            const clickable = canSelect && !isFixed
            const rowClasses = cn(
              isFixed && 'bg-slate-50',
              isSelected && 'bg-sky-100',
              clickable && 'cursor-pointer hover:bg-sky-200/40',
            )
            const fixedCoordClass = isFixed
              ? 'bg-slate-50 text-slate-500'
              : 'text-slate-600'
            return (
              <TableRow
                key={`${row.symbol}-${row.index}`}
                className={rowClasses}
                onClick={
                  clickable && onRowClick
                    ? () => onRowClick(row.index)
                    : undefined
                }
              >
                {showIndex ? (
                  <TableCell className="text-xs font-medium text-slate-500">
                    {row.index + 1}
                  </TableCell>
                ) : null}
                <TableCell className="font-semibold text-slate-700">
                  {row.symbol}
                </TableCell>
                <TableCell className={cn('font-mono', fixedCoordClass)}>
                  {formatNumber(row.x, digits)}
                </TableCell>
                <TableCell className={cn('font-mono', fixedCoordClass)}>
                  {formatNumber(row.y, digits)}
                </TableCell>
                <TableCell className={cn('font-mono', fixedCoordClass)}>
                  {formatNumber(row.z, digits)}
                </TableCell>
              </TableRow>
            )
          })
        ) : (
          <TableRow>
            <TableCell
              colSpan={colSpan}
              className="py-10 text-center text-sm text-slate-500"
            >
              {emptyText}
            </TableCell>
          </TableRow>
        )}
      </TableBody>
    </Table>
  )
}
